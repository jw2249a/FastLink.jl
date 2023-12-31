function pool_lookup_table(vec::PooledVector)::Dict{UInt32, Vector{UInt32}}
    lookup_by_id = Dict{UInt32, Vector{UInt32}}()
    for name in collect(vec.invpool.vals)
        lookup_by_id[name] = findall(x -> x == name, vec.refs) 
    end
    return lookup_by_id
end

mutable struct CandidateLetterInfo
    name_index::UInt32
    len::UInt8
    mask::UInt16
end

mutable struct CandidateScore
    matches::UInt8
    used::UInt16
    used_exact::UInt16
    len_partial::UInt16
    last_match_letter_index::UInt16
    transposition_count::UInt8

    function CandidateScore(len::Number)
        new(0, 0, 0, UInt16(1024.0 ÷ len), 0, 0)
    end
end



function build_candidate_lookup(name_vec::Vector{T};
                                allcaps::Bool = true)::Vector{Vector{CandidateLetterInfo}} where T <: AbstractString
    # defining the range of all the letters
    firstletter::UInt8 = allcaps ? 0x41 : 0x61
    letters::UnitRange{UInt8} = UnitRange(firstletter,firstletter+0x19)
    
    letter_lookup = Vector{Vector{CandidateLetterInfo}}(undef, 26)
    
    Threads.@threads for letter in letters
        candidate_infos = Vector{CandidateLetterInfo}()
        for (name_index, name) in enumerate(name_vec)
                mask = UInt16(0)
                for (matching_index_in_name, c) in enumerate(codeunits(name))
                    if c == letter
                        mask += UInt16(2)^(matching_index_in_name - 1)
                    end
                end
                if mask != 0
                    push!(candidate_infos, CandidateLetterInfo(name_index, UInt8(ncodeunits(name)), mask))
                end

            letter_lookup[letter - firstletter + 1] = candidate_infos
        end
    end
    return letter_lookup
end

function build_candidate_lookup(name_vec::Vector{Union{Missing,T}};
                                allcaps::Bool = true)::Vector{Vector{CandidateLetterInfo}} where T <: AbstractString
    # defining the range of all the letters
    firstletter::UInt8 = allcaps ? 0x41 : 0x61
    letters::UnitRange{UInt8} = UnitRange(firstletter,firstletter+0x19)
    
    letter_lookup = Vector{Vector{CandidateLetterInfo}}(undef, 26)
    
    Threads.@threads for letter in letters
        candidate_infos = Vector{CandidateLetterInfo}()
        for (name_index, name) in enumerate(name_vec)
            if !ismissing(name)
                mask = UInt16(0)
                for (matching_index_in_name, c) in enumerate(codeunits(name))
                    if c == letter
                        mask += UInt16(2)^(matching_index_in_name - 1)
                    end
                end
                if mask != 0
                    push!(candidate_infos, CandidateLetterInfo(name_index, UInt8(ncodeunits(name)), mask))
                end
            end
            letter_lookup[letter - firstletter + 1] = candidate_infos
        end
    end
    return letter_lookup
end


function get_query_mask(cl::Int,query_mask::UInt16, min_match_dist::Int)::UInt16
    match_distance=cl <= 3 || cl ÷ 2 - 1 <= min_match_dist ? min_match_dist : cl ÷ 2 - 1
    for _ in 0:match_distance
        query_mask = query_mask << 1 | query_mask;
        query_mask = query_mask >> 1 | query_mask;
    end
    return query_mask;
end
"""Converts a string to a list of bitmasks"""
function maskify(query::T;
                 space_char::UInt8=0x60)::Vector{Tuple{UInt8, Vector{UInt16}}} where T <: AbstractString
    min_match_dist = let len=ncodeunits(query); len > 3 ? len ÷ 2 - 1 : 0 end

    [(c-space_char,
      [get_query_mask(cl,0x0001 << i,min_match_dist) for cl in 1:16])
    for (i, c) in enumerate(replace(codeunits(query), 0x20 => space_char))]
end

function score_letter!(candidate_score::CandidateScore, query_mask::UInt16, candidate_mask::UInt16, query_index::Int)
    whole_mask_result = query_mask & candidate_mask
    check_used_result = xor(whole_mask_result | candidate_score.used, candidate_score.used)
    last_match_letter_index = 0x0001 << trailing_zeros(check_used_result)
    mask_result = check_used_result & last_match_letter_index
    is_match_mask = UInt16(trailing_ones(mask_result >> trailing_zeros(mask_result)))
    candidate_score.used |= mask_result
    candidate_score.used_exact |= mask_result & (0x0001 << query_index-1)
    candidate_score.matches += is_match_mask
    candidate_score.transposition_count += UInt8(mask_result - 1 < candidate_score.last_match_letter_index)
    candidate_score.last_match_letter_index |= mask_result
end

function calculate_jaro_winkler(score::CandidateScore, query_partial::UInt16; p=0.1)::Float64
    transpositions = score.transposition_count > score.matches / 2 ? score.transposition_count - 1 : score.transposition_count
    partial = ((score.matches * score.len_partial) + (score.matches * query_partial)) / 1024.0
    jaro = (partial + 1.0 - (transpositions / score.matches)) / 3.0
    l = trailing_ones(UInt16(score.used_exact & 0x000f))
    return jaro + p * l * (1.0 - jaro)
end

function format_vector(x)
    first(replace(replace(x, " "=>""), !isascii=>""),16)
end

function gammaCKfuzzy!(vecA::PooledVector,vecB::PooledVector, results::SubArray, array_2Dindex::Function, dims::Tuple;
                       cut_a=0.92, cut_b=0.88,upper=true)
    match2 = [true, true]
    match1 = [true, false]
    missingval = [false, true]

    vecA = Vector(map(x-> ismissing(x) || format_vector(x)=="" ? missing : format_vector(x), vecA))
    vecB = Vector(map(x-> ismissing(x) || format_vector(x)=="" ? missing : format_vector(x), vecB))
    
    missingvals_x=filter(x->ismissing(x),vecA)
    missingvals_y=filter(x->ismissing(x),vecB)

    filter!(x->!ismissing(x),vecA)
    filter!(x->!ismissing(x),vecB)
    
    lookup_a_by_name = Dict(name => findall(x -> x == name, vecA) for name in vecA)
    copy!(vecA,sort!(unique(vecA)))
    lookup_a_by_id = Dict(i => lookup_a_by_name[name] for (i, name) in enumerate(vecA))

    lookup_b_by_name = Dict(name => findall(x -> x == name, vecB) for name in vecB)
    copy!(vecB,sort!(unique(vecB)))
    lookup_b_by_id = Dict(i => lookup_b_by_name[name] for (i, name) in enumerate(vecB))
    
    base_candidate_lookup = build_candidate_lookup(vecB)
    base_candidate_scores = [CandidateScore(ncodeunits(name)) for name in vecB]
    
    Threads.@threads for (new_a_id,query_name) in collect(enumerate(vecA))
        query_masks_lookup = maskify(query_name,space_char=0x40)
        query_partial = UInt16(1024 ÷ ncodeunits(query_name))
        candidate_scores = deepcopy(base_candidate_scores)
        for (query_index, (letter_index, query_mask_by_candidate_len)) in enumerate(query_masks_lookup)
            for c_info in @inbounds base_candidate_lookup[letter_index]
                candidate_score = candidate_scores[c_info.name_index]
                query_mask = query_mask_by_candidate_len[c_info.len]
                score_letter!(candidate_score, query_mask, c_info.mask, query_index)
            end
        end

        a_ids = lookup_a_by_id[new_a_id]
        for (score_i, score) in enumerate(candidate_scores)
            # if present calculate scores
            jw = calculate_jaro_winkler(score, query_partial)
            if jw >= cut_a
                b_ids = lookup_b_by_id[score_i]
                for a_id in a_ids
                    for b_id in b_ids
                        results[array_2Dindex(Int(a_id),Int(b_id)),:] = match2
                    end
                end
            elseif jw >= cut_b
                b_ids = lookup_b_by_id[score_i]
                for a_id in a_ids
                    for b_id in b_ids
                        results[array_2Dindex(Int(a_id),Int(b_id)),:] = match1
                    end
                end
            end
            
        end

    end

    if length(missingvals_x) != 0
        missingindices = findall(ismissing.(vecA))
        Threads.@threads for iy in 1:dims[2]
            for ix in missingindices
                results[array_2Dindex(ix,iy),:] = missingval
            end
        end
    end
    # set all to missing where y is missing
    if length(missingvals_y) != 0
        missingindices = findall(ismissing.(vecB))
        Threads.@threads for ix in 1:dims[1]
            for iy in missingindices
                results[array_2Dindex(ix,iy),:] = missingval
            end
        end
    end
    
    return nothing
end

function pool_lookup_table(vec::Vector{UInt32},len::UInt32)
    lookup_by_id = fill(Vector{UInt32}(), len)
    Threads.@threads for pool_index in UInt32(1):len
        lookup_by_id[pool_index] = (1:length(vec))[vec .=== pool_index]
    end
    return lookup_by_id
end

mutable struct CandidateLetterInfo
    name_index::UInt32
    len::UInt8    
    mask::UInt16
end
"""struct to hold candidate comparisons"""
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

function build_candidate_scores(x::Vector)
    [CandidateScore(ismissing(name) ? 1024 : min(ncodeunits(name),16)) for name in x]
end

function build_candidate_lookup(name_vec::Vector;
                                spaceletter::UInt8 = 0x40,lastletter::UInt8 = 0x5a)::Vector{Vector{CandidateLetterInfo}}
    # defining the range of all the letters
    letters::UnitRange{UInt8} = UnitRange(spaceletter+1,lastletter)

    letter_lookup = fill(Vector{CandidateLetterInfo}(),lastletter-spaceletter+1)
    #letter_lookup = Vector{Vector{CandidateLetterInfo}}(undef, lastletter-spaceletter+1)
    
    for letter in letters
        candidate_infos = Vector{CandidateLetterInfo}()
        for (name_index, name) in collect(enumerate(name_vec))
            if !ismissing(name)
                mask = UInt16(0)
                name=first(name,16)
                for (matching_index_in_name, c) in enumerate(codeunits(name))
                    if c == letter
                        mask += UInt16(2)^(matching_index_in_name - 1)
                    end
                end
                if mask != 0
                    push!(candidate_infos, CandidateLetterInfo(name_index, UInt8(ncodeunits(name)), mask))
                end
            end
            letter_lookup[letter - spaceletter + 1] = candidate_infos
        end
    end
    # alter it for spaces
    letter=0x20
    candidate_infos = Vector{CandidateLetterInfo}()
    for (name_index, name) in collect(enumerate(name_vec))
        if !ismissing(name)
            mask = UInt16(0)
            name=first(name,16)
            for (matching_index_in_name, c) in enumerate(codeunits(name))
                if c == letter
                    mask += UInt16(2)^(matching_index_in_name - 1)
                end
            end
            if mask != 0
                push!(candidate_infos, CandidateLetterInfo(name_index, UInt8(ncodeunits(name)), mask))
            end
        end
        letter_lookup[1] = candidate_infos
    end
    return letter_lookup
end

function get_query_mask(cl::UInt8,query_mask::UInt16, min_match_dist::UInt8)::UInt16
    match_distance=cl <= 3 || cl ÷ 2 - 1 <= min_match_dist ? min_match_dist : cl ÷ 2 - 1
    for _ in 0:match_distance
        query_mask = query_mask << 1 | query_mask;
        query_mask = query_mask >> 1 | query_mask;
    end
    return query_mask;
end
"""Converts a string to a list of bitmasks"""
function maskify(query::T,len::UInt8;
                 space_char::UInt8=0x40, max_char::UInt8=0x5a)::Vector{Tuple{UInt8, Vector{UInt16}}} where T <: AbstractString
    min_match_dist = UInt8(len > 3 ? len ÷ 2 - 1 : 0)
    
    [(c-space_char+UInt8(1),
      [get_query_mask(cl,0x0001 << i,min_match_dist) for cl in UInt8(1):UInt8(16)])
    for (i, c) in enumerate(filter(x-> x>=space_char && x<=max_char, replace(codeunits(query), 0x20 => space_char)))]
end
"""scores individual letters in based on two candidate scores"""
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

function find_missing_index(d::Dict)::UInt32
    return let x = [i for i in values(filter(x->ismissing(x[1]), d))]; length(x) > 0 ? x[1] : UInt32(0) ; end
end

"""
Fuzzy string comparison of two columns (with 2 levels of similarity).

Fuzzy JW based on:
https://tech.popdata.org/speeding-up-Jaro-Winkler-with-rust-and-bitwise-operations/

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::SubArray`: ResultMatrix object's result_matrix.
- `array_2Dindex::Function`: ResultMatrix object's array_2Dindex function
- `dims::Tuple`: ResultMatrix object's dims.
- `cut_a::Float=0.92`: Lower bound for close string distances.
- `cut_b::Float=0.88`: Lower bound for partial string distance.
- `upper::Bool=true`: Whether input string is uppercase.
- `w`: Winkler weight for jw string distance.
"""
function gammaCKfuzzy!(vecA::PooledVector,vecB::PooledVector, results::SubArray, array_2Dindex::Function, dims::Tuple; cut_a=0.92, cut_b=0.88,upper=true,w=0.1)

    # functions that update the results view
    function update_results!(a_ids::Vector{UInt32},
                             b_ids::Vector{UInt32},
                             val::Matrix{Bool};
                             results::SubArray=results,
                             array_2Dindex::Function=array_2Dindex)
        results[[array_2Dindex(UInt32(ia),UInt32(ib)) for ia in a_ids for ib in b_ids],:] .= val
        return nothing
    end
    function update_results!(a_ids::Vector{UInt32},
                             b_ids::UnitRange{UInt32},
                             val::Matrix{Bool};
                             results::SubArray=results,
                             array_2Dindex::Function=array_2Dindex)
        results[[array_2Dindex(UInt32(ia),UInt32(ib)) for ia in a_ids for ib in b_ids],:] .= val
        return nothing
    end

    
    match2 = [true  true]
    match1 = [true false]
    missingval = [false true]

    if upper
        space_char,max_char = 0x40,0x5a
    else
        space_char,max_char = 0x60,0x7a
    end
    
    lenA = UInt32(length(vecA.pool))
    lenB = UInt32(length(vecB.pool))

    lookup_a_by_id=pool_lookup_table(vecA.refs, lenA)
    lookup_b_by_id=pool_lookup_table(vecB.refs, lenB)
    
    missingindexA = find_missing_index(vecA.invpool)
    
    base_candidate_lookup = build_candidate_lookup(vecB.pool,spaceletter=space_char,lastletter=max_char)
    base_candidate_scores = build_candidate_scores(vecB.pool)
    
    Threads.@threads for (query_name,new_a_id) in collect(vecA.invpool)
        # pass if query is missing val
        if new_a_id === missingindexA
            update_results!(lookup_a_by_id[new_a_id],UInt32(1):UInt32(dims[2]),missingval)
            continue
        end
        
        query_len = UInt8(min(ncodeunits(query_name),16))
        query_masks_lookup = maskify(query_name,query_len,space_char=0x40,max_char=max_char)
        query_partial = UInt16(1024 ÷ query_len)
        candidate_scores = deepcopy(base_candidate_scores)
        
        for (query_index, (letter_index, query_mask_by_candidate_len)) in enumerate(query_masks_lookup)
            for c_info in base_candidate_lookup[letter_index]
                candidate_score = candidate_scores[c_info.name_index]
                query_mask = query_mask_by_candidate_len[c_info.len]
                score_letter!(candidate_score, query_mask, c_info.mask, query_index)
            end
        end
        
        a_ids = lookup_a_by_id[new_a_id]
        for (score_i, score) in enumerate(candidate_scores)
            if score.len_partial === UInt16(1)                
                update_results!(a_ids,lookup_b_by_id[score_i],missingval)
                continue
            end
            # if present calculate scores
            jw = calculate_jaro_winkler(score, query_partial,p=w)
            if jw >= cut_a
                update_results!(a_ids,lookup_b_by_id[score_i],match2)
            elseif jw >= cut_b
                update_results!(a_ids,lookup_b_by_id[score_i],match1)
            end         
        end
    end

    return nothing
end


"""
Fuzzy string comparison of two columns (with 1 level of similarity).

Fuzzy JW based on:
https://tech.popdata.org/speeding-up-Jaro-Winkler-with-rust-and-bitwise-operations/

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::SubArray`: ResultMatrix object's result_matrix.
- `array_2Dindex::Function`: ResultMatrix object's array_2Dindex function
- `dims::Tuple`: ResultMatrix object's dims.
- `cut_a::Float=0.92`: String distance cutoff.
- `upper::Bool=true`: Whether input string is uppercase.
- `w`: Winkler weight for jw string distance.
"""
function gammaCK2fuzzy!(vecA::PooledVector,vecB::PooledVector, results::SubArray, array_2Dindex::Function, dims::Tuple; cut_a=0.92, upper=true,w=0.1)

    # functions that update the results view
    function update_results!(a_ids::Vector{UInt32},
                             b_ids::Vector{UInt32},
                             val::Matrix{Bool};
                             results::SubArray=results,
                             array_2Dindex::Function=array_2Dindex)
        results[[array_2Dindex(UInt32(ia),UInt32(ib)) for ia in a_ids for ib in b_ids],:] .= val
        return nothing
    end
    function update_results!(a_ids::Vector{UInt32},
                             b_ids::UnitRange{UInt32},
                             val::Matrix{Bool};
                             results::SubArray=results,
                             array_2Dindex::Function=array_2Dindex)
        results[[array_2Dindex(UInt32(ia),UInt32(ib)) for ia in a_ids for ib in b_ids],:] .= val
        return nothing
    end

    
    match2 = [true true]
    missingval = [false true]

    if upper
        space_char,max_char = 0x40,0x5a
    else
        space_char,max_char = 0x60,0x7a
    end
    
    lenA = UInt32(length(vecA.pool))
    lenB = UInt32(length(vecB.pool))

    lookup_a_by_id=pool_lookup_table(vecA.refs, lenA)
    lookup_b_by_id=pool_lookup_table(vecB.refs, lenB)
    
    missingindexA = find_missing_index(vecA.invpool)
    
    base_candidate_lookup = build_candidate_lookup(vecB.pool,spaceletter=space_char,lastletter=max_char)
    base_candidate_scores = build_candidate_scores(vecB.pool)
    
    Threads.@threads for (query_name,new_a_id) in collect(vecA.invpool)
        # pass if query is missing val
        if new_a_id === missingindexA
            update_results!(lookup_a_by_id[new_a_id],UInt32(1):UInt32(dims[2]),missingval)
            continue
        end
        
        query_len = UInt8(min(ncodeunits(query_name),16))
        query_masks_lookup = maskify(query_name,query_len,space_char=0x40,max_char=max_char)
        query_partial = UInt16(1024 ÷ query_len)
        candidate_scores = deepcopy(base_candidate_scores)
        
        for (query_index, (letter_index, query_mask_by_candidate_len)) in enumerate(query_masks_lookup)
            for c_info in base_candidate_lookup[letter_index]
                candidate_score = candidate_scores[c_info.name_index]
                query_mask = query_mask_by_candidate_len[c_info.len]
                score_letter!(candidate_score, query_mask, c_info.mask, query_index)
            end
        end
        
        a_ids = lookup_a_by_id[new_a_id]
        for (score_i, score) in enumerate(candidate_scores)
            if score.len_partial === UInt16(1)                
                update_results!(a_ids,lookup_b_by_id[score_i],missingval)
                continue
            end
            # if present calculate scores
            jw = calculate_jaro_winkler(score, query_partial,p=w)
            if jw >= cut_a
                update_results!(a_ids,lookup_b_by_id[score_i],match2)
            end         
        end
    end

    return nothing
end

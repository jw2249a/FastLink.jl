module patterns

using ..matchpatterns, ..emlink, ..tf, ..DiBitMat
import ..match2, ..nonmatch

export indices_to_uids, process_comparisons
function indices_to_uids(vecA, vecB,
                         indices::Vector{Vector{ComparisonIndex}}
                         )
    batch_size=500
    inds=eachindex(indices)
    paired_ids = [Vector{Tuple}() for _ in inds]
    Threads.@threads for i in inds
        len=length(indices[i])
        lk = ReentrantLock()
        Threads.@threads for first_val in 1:batch_size:len
            local_paired_ids=Vector{Tuple}()
            last_val = min(first_val + batch_size - 1, len)
            for ii in first_val:last_val
                push!(local_paired_ids,(vecA[indices[i][ii].row],vecB[indices[i][ii].col]))
            end
            lock(lk) do
                append!(paired_ids[i],local_paired_ids)
            end
        end
    end
    return paired_ids
end

# Decomposes match vector index back into original indices from initial inputs.
function get_2Dindex(index::T, nrows::Int) where T <: Integer
    zero_based_index = index - 1
    row = Int(mod(zero_based_index, nrows)) + 1
    col = Int(div(zero_based_index, nrows)) + 1
    return ComparisonIndex(row, col)
end

function get_local_patterns(x::Vector{Vector{UInt8}}, N::Int, S::Int)
    patterns=Vector{Vector{UInt8}}()
    hashes=Vector{UInt64}()
    indices=Vector{Vector{UInt16}}()

    for i in 1:S
        pattern=zeros(UInt8,N)
        for n in 1:N
            pattern[n]=x[n][i]
        end
        pattern_hash=hash(pattern)
        id = findfirst(pattern_hash .=== hashes)
        if isnothing(id)
            push!(patterns,pattern)
            push!(hashes,pattern_hash)
            push!(indices,[i])
        else
            push!(indices[id],i)
        end
    end
    return LocalPatterns(patterns,indices,hashes)
end


function get_match_patterns(res::Vector{DiBitMatrix})
    matches=MatchPatterns()
    N = length(res)
    dimy=res[1].nrows
    len=Int(res[1].data.len)
    lk = ReentrantLock()
    Threads.@threads for first_loc in 0:1024:len
        last_loc = first_loc + 1024
        if last_loc > len
            last_loc=len
        end
        x=[res[n].data[(first_loc+1):last_loc] for n in 1:N]
        patterns=get_local_patterns(x,N,last_loc-first_loc)
        for i in eachindex(patterns.hashes)
            lock(lk) do 
                id = findfirst(patterns.hashes[i] .=== matches.hashes)
                if isnothing(id)
                    push!(matches.patterns,patterns.patterns[i])
                    push!(matches.hashes,patterns.hashes[i])
                    push!(matches.indices,get_2Dindex.(first_loc .+ patterns.indices[i],dimy))
                else
                    append!(matches.indices[id],get_2Dindex.(first_loc .+ patterns.indices[i],dimy))
                end
            end
        end
    end
    return matches
end

function get_match_patterns(res::Vector{DiBitMatrix}, tf_tables::Dict{String, Vector{Vector{Float16}}},
                            tf_vars::Vector{String}, tf_indices::Vector{Int64}, isexact=Bool[])
    tf_patterns = Dict("relevant_tf_indices"=>Vector{Int64}[],
                       "tf_denom_vals"=>Vector{Vector{Float16}}[])  
    matches=MatchPatterns()
    N = length(res)
    dimy=res[1].nrows
    len=Int(res[1].data.len)
    lk = ReentrantLock()
    Threads.@threads for first_loc in 0:1024:len
        last_loc = first_loc + 1024
        if last_loc > len
            last_loc=len
        end
        x=[res[n].data[(first_loc+1):last_loc] for n in 1:N]
        patterns=get_local_patterns(x,N,last_loc-first_loc)
        for i in eachindex(patterns.hashes)
            lock(lk) do
                id = findfirst(patterns.hashes[i] .=== matches.hashes)
                pattern_indices = get_2Dindex.(first_loc .+ patterns.indices[i],dimy)
                if isnothing(id)
                    push!(matches.patterns,patterns.patterns[i])
                    push!(matches.hashes,patterns.hashes[i])
                    push!(matches.indices, pattern_indices)
                    relevant_tf_indices = find_tf_pattern_vars(patterns.patterns[i], tf_indices)
                    tfi_ids = [findfirst(tfi .== tf_indices) for tfi in relevant_tf_indices]
                    push!(tf_patterns["relevant_tf_indices"], relevant_tf_indices)
                    push!(tf_patterns["tf_denom_vals"], [match_level_tf_lookup(tf_tables[tf_vars[tfi]], pattern_indices, isexact[tfi]) for tfi in tfi_ids])
                else
                    tfi_ids = [findfirst(tfi .== tf_indices) for tfi in tf_patterns["relevant_tf_indices"][id]]
                    append!(matches.indices[id], pattern_indices)
                    for (tfi_loc, tfi) in enumerate(tfi_ids)
                        append!(tf_patterns["tf_denom_vals"][id][tfi_loc], match_level_tf_lookup(tf_tables[tf_vars[tfi]], pattern_indices, isexact[tfi]))
                    end
                end
            end
        end
    end
    return (tf_patterns, matches)
end



function match_level_tf_lookup(tf_table::Vector{Vector{Float16}}, pattern_indices::Vector{ComparisonIndex}, isexact::Bool)
    if isexact
        calculate_tf_denom = calculate_tf_denom_exact
    else
        calculate_tf_denom = calculate_tf_denom_fuzzy
    end
    return [calculate_tf_denom(tf_table[1][pi.row],tf_table[2][pi.col]) for pi in pattern_indices]
end

function get_match_patternids(resultsEM; base="log2")
    if base == "log2"
        bf_to_prob = bf2_to_probability
    else
        bf_to_prob = bf_to_probability
    end

    resultsEM["patterns_w"].weights .|> Float64 .|>
        bf_to_prob .|> (x -> x >= resultsEM["threshold_match"]) |>
        findall
end

function update_pattern_level_DiBit!(pattern::DiBitMatrix,
                                     match_vector::BitVector,
                                     indices::Vector{ComparisonIndex})
    count = length(match_vector)
    Threads.@threads for i in collect(1:count)
        pattern[indices[i]] = match_vector[i] ? match2 : nonmatch
    end
    return nothing
end

function patterns_to_DiBit(resultsTF::Vector{Dict{String, Any}},
                           comparisons::Vector{Vector{ComparisonIndex}},
                           _dims::Tuple{Int64,Int64})
    patterns = DiBitMatrix(_dims...)

    # find non adjusted matches
    pattern_ids = findall([any(p["ismatch"]) && p["tf_adjusted"] == false for p in resultsTF])
    Threads.@threads for pattern_id in pattern_ids
        Threads.@threads for obs in comparisons[pattern_id]
            patterns[obs] = match2
        end
    end

    # find tf_adjusted matches
    pattern_ids = findall([any(p["ismatch"]) && p["tf_adjusted"] for p in resultsTF])
    Threads.@threads for pattern_id in pattern_ids
        update_pattern_level_DiBit!(patterns,
                                    resultsTF[pattern_id]["ismatch"],
                                    comparisons[pattern_id])
    end

    return patterns
end

function match_and_link(patterns::Vector{DiBitMatrix}, e::Dict{String, Any}, _dims::Tuple{Int64,Int64},
                        final_name::String)
    @info "getting match patterns"
    counts=get_match_patterns(patterns)
    @info "running emlink"
    resultsEM=emlinkMARmov(counts.patterns,
                           length.(counts.indices),
                           _dims,
                           e["variables"];
                           e["parameters"]...)
    if e["name"] != final_name
        pattern_ids = get_match_patternids(resultsEM; base="log")
        resultsTF = [Dict{String,Any}("ismatch" => pi ∈ pattern_ids,"tf_adjusted" => false)
                     for pi in 1:resultsEM["number_of_unique_patterns"]]
        return patterns_to_DiBit(resultsTF, counts.indices, _dims)
    else
        return (counts, resultsEM)
    end
end


function match_and_link(patterns::Vector{DiBitMatrix}, e::Dict{String, Any}, _dims::Tuple{Int64,Int64},
                        final_name::String, tf_vars::Vector{String}, parameters::Dict{String, Dict{String, Any}},
                        tf_tables::Dict{String, Vector{Vector{Float16}}})
    tf_indices=[findfirst(i .== e["variables"]) for i in tf_vars]
    isexact = [parameters[varname]["method"] == "exact" for varname in tf_vars]
    tf_patterns, counts = get_match_patterns(patterns,tf_tables,tf_vars, tf_indices, isexact)
    resultsEM = emlinkMARmov(counts.patterns,
                             length.(counts.indices),
                             _dims,
                             e["variables"];
                             e["parameters"]...)
    
    tf_prior_weights = get_tf_adjustment_prior_weights(parameters, tf_vars)
    resultsTF = generate_tf_adjustment_dict(resultsEM, e, tf_vars, tf_patterns, tf_prior_weights; base="log")

    if e["name"] != final_name
        return patterns_to_DiBit(resultsTF, counts.indices, _dims)
    else
        return (counts, resultsEM, resultsTF)
    end
end

function process_comparisons(res::Dict{String, DiBitMatrix},
                             emlink_configuration::Vector{Vector{Dict{String, Any}}},
                             _dims::Tuple{Int64,Int64},
                             parameters::Dict{String, Dict{String, Any}},
                             tf_tables::Dict{String, Vector{Vector{Float16}}})
    final_name = last(emlink_configuration)[1]["name"]
    for emconfig in emlink_configuration
        for e in emconfig
            tf_vars=intersect(e["variables"],keys(tf_tables))
            patterns=[pop!(res, v) for v in e["variables"]]
            if isempty(tf_vars)
                if e["name"] != final_name
                    res[e["name"]] = match_and_link(patterns, e, _dims, final_name)
                else
                    return match_and_link(patterns, e, _dims, final_name)
                end
                
            else
                if e["name"] != final_name
                    res[e["name"]] = match_and_link(patterns, e, _dims, final_name, tf_vars, parameters, tf_tables)
                else
                    return match_and_link(patterns, e, _dims, final_name, tf_vars, parameters, tf_tables)
                end
                
            end
        end
    end
end
end #module patterns

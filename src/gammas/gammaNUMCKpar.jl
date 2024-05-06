function allow_missing(x::Vector)
    return copy!(Vector{Union{eltype(x), Missing}}(),x)
end

# difference categorization functions functions 
function get_diff1(x::T,y::T,cut_a::T)::UInt8 where T <: Number
    y - x < cut_a ? match2 : nonmatch
end
function get_diff2(x::T,y::T,cut_a::T, cut_b::T)::UInt8 where T <: Number
    d = y - x
    if d <= cut_a
        return match2
    else
        if d <= cut_b
            return match1
        else
            return nonmatch
        end
    end
end

function get_diff2(x::T, y::Missing,cut_a::T, cut_b::T)::UInt8 where T <: Number
    return missingval
end


fix_id(id::Int, N_a::Int) = id <= N_a ? (id, true) : (id - N_a, false)


"""
Numeric comparison of two columns

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `res::DiBitMatrix`: DiBitMatrix.
- `cut_a::Number=1`: Lower bound for close string distances.
- `cut_b::Number=2`: Lower bound for partial string distances.
"""
function gammaNUMCKpar!(vecA::Vector, vecB::Vector,
                        results::DiBitMatrix;
                        cut_a=1,cut_b=2,
                        partial::Bool=true)
    @info "running gammanumckpar"
    N_a = length(vecA)
    N_b = length(vecB)

    # whether it is two levels or 1
    if partial
        get_diff = get_diff2
    else
        get_diff = get_diff1
    end

    vecA=allow_missing(vecA)
    vecB=allow_missing(vecB)

    # coercing cuts in case they are wrong type
    if (eltype(vecA) <: Union{Integer,Missing}) & (eltype(vecB) <: Union{Integer,Missing})
        cut_a = Int(cut_a)
        cut_b = Int(cut_b)
    else
        cut_a = Float64(cut_a)
        cut_b = Float64(cut_b)
    end

    # get the sorted indices of a large copied array
    append!(vecA,vecB)
    vecC=sortperm(vecA)
    copy!(vecA,vecA[vecC])
    len=length(vecA)
    
    # preallocation of ranges for the threads
    tids = Threads.nthreads()
    breaksize= len ÷ (tids-1)
    starts = (collect(0:(tids-1))) .* breaksize .+ 1
    ends = starts[:]
    if last(starts) == len
        pop!(starts)
        popfirst!(ends)
    else
        ends = append!(starts[:], [len])
        popfirst!(ends)
    end
    tids=length(starts)
    
    Threads.@threads for tid in 1:tids    
        start=starts[tid]
        s_end=ends[tid]

        under_consideration_value = [vecA[start]]
        current_id, tableA_bool = fix_id(vecC[start], N_a)
        under_consideration_id = [current_id]
        tableA_bools = [tableA_bool]
        match_values = Vector{UInt8}()
        
        for i in 1:(s_end-start)
            itarg=start+i
            obs=vecA[itarg]
            current_id, tableA_bool = fix_id(vecC[itarg], N_a)
            match_values=get_diff.(under_consideration_value, obs, cut_a, cut_b)
            if last(match_values) === missingval
                for ii in itarg:s_end
                    current_id, tableA_bool = fix_id(vecC[ii], N_a)
                    if tableA_bool
                        for icol in 1:N_b
                            results[current_id, icol] = missingval
                        end
                    else
                        for irow in 1:N_a
                            results[irow, current_id] = missingval
                        end
                    end
                end
                under_consideration_id=Int64[]
                under_consideration_value=Int64[]
                tableA_bools = Bool[]
                match_values = Vector{UInt8}()
                break # break if we run into missing values
            end
            removed_values=0
            for m in eachindex(match_values)
                if match_values[m] === nonmatch
                    popfirst!(under_consideration_id)
                    popfirst!(under_consideration_value)
                    popfirst!(tableA_bools)
                    removed_values += 1
                elseif match_values[m] === match2
                    # iterate through and break because they all match
                    for idc in (m-removed_values):(length(match_values)-removed_values) 
                        if tableA_bool ⊻ tableA_bools[idc]
                            if tableA_bool
                                results[current_id,under_consideration_id[idc]] = match2
                            else
                                results[under_consideration_id[idc], current_id] = match2
                            end
                        end
                    end
                    break
                else
                    if tableA_bool ⊻ tableA_bools[m-removed_values]
                        if tableA_bool
                            results[current_id,under_consideration_id[m-removed_values]] = match1
                        else
                            results[under_consideration_id[m-removed_values], current_id] = match1
                        end
                    end
                end
            end
            push!(under_consideration_value, obs)
            push!(under_consideration_id, current_id)
            push!(tableA_bools,tableA_bool)
        end

        i = (s_end-start)
        
        # clear out the overlap
        while length(under_consideration_id) > 0
            i += 1
            itarg=start+i
            if itarg === len
                break
            end
            obs=vecA[itarg]
            current_id, tableA_bool = fix_id(vecC[itarg], N_a)
            match_values=get_diff.(under_consideration_value, obs, cut_a, cut_b)

            removed_values=0
            for m in eachindex(match_values)
                if match_values[m] === nonmatch
                    popfirst!(under_consideration_id)
                    popfirst!(under_consideration_value)
                    popfirst!(tableA_bools)
                    removed_values += 1
                elseif match_values[m] === match2
                    for idc in (m-removed_values):(length(match_values)-removed_values) # iterate through and break because they all match
                        if tableA_bool ⊻ tableA_bools[idc]
                            if tableA_bool
                                results[current_id,under_consideration_id[idc]] = match2
                            else
                                results[under_consideration_id[idc], current_id] = match2
                            end
                        end
                    end
                    break
                else
                    if tableA_bool ⊻ tableA_bools[m-removed_values]
                        if tableA_bool
                            results[current_id,under_consideration_id[m-removed_values]] = match1
                        else
                            results[under_consideration_id[m-removed_values], current_id] = match1
                        end
                    end
                end
            end
        end
        
    end
    @info "finish gammanumckpar"
    return nothing
end




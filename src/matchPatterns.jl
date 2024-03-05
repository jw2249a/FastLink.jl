struct ComparisonIndex
    row::UInt32
    col::UInt32
    function ComparisonIndex(row::Int,col::Int)
        new(UInt32(row),UInt32(col))
    end
end

struct LocalPatterns
    patterns::Vector{Vector{UInt8}}
    indices::Vector{Vector{Int}}
    hashes::Vector{UInt64}
end

struct MatchPatterns
    patterns::Vector{Vector{UInt8}}
    indices::Vector{Vector{ComparisonIndex}}
    hashes::Vector{UInt64}

    function MatchPatterns()
        new(Vector{Vector{UInt8}}(), Vector{ComparisonIndex}(),Vector{UInt64}())
    end
end

function indices_to_uids(vecA, vecB,
                                    indices::Vector{Vector{ComparisonIndex}}
                                 ) 
    inds=eachindex(indices)
    paired_ids = [Vector{Tuple}() for _ in inds]
    Threads.@threads for i in inds
        len=length(indices[i])
        lk = ReentrantLock()
        Threads.@threads for first_val in 1:500:len
            local_paired_ids=Vector{Tuple}()
            last_val = first_val + min(len-first_val,500-first_val)
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

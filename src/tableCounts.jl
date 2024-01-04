function generate_bitvector(tid::Int, vlen::Int)::Vector{Bool}
    return [(((tid-1) >> j) & 1) == 1 for j in 0:(vlen - 1)]
end

function find_matches(res_matrix::SubArray,id::Int,vlen::Int)
    targ_vec=generate_bitvector(id, vlen)
    match_count::UInt64=0
    row_count::UInt64=0
    matches::Vector{UInt64}=[]
    for i in eachrow(res_matrix)
        row_count += 1
        if i == targ_vec
            match_count += 1
            append!(matches, row_count)
        end
    end
    return (targ_vec, match_count, matches)
end

function find_matches(res_matrix::Matrix,id::Int,vlen::Int)
    targ_vec=generate_bitvector(id, vlen)
    match_count::UInt64=0
    row_count::UInt64=0
    matches::Vector{UInt64}=[]
    for i in eachrow(res_matrix)
        row_count += 1
        if i == targ_vec
            match_count += 1
            append!(matches, row_count)
        end
    end
    return (targ_vec, match_count, matches)
end

function tableCounts_bf(res_matrix::SubArray)
    vlen=size(res_matrix)[2]
    nlevels=2^vlen
    out=fill((falses(0),UInt64(0),zeros(UInt64,0)),nlevels)
    Threads.@threads for i in 1:nlevels
        out[i]=find_matches(res_matrix,i, vlen)
    end
    return out
end

function tableCounts_bf(res_matrix::Matrix)
    vlen=size(res_matrix)[2]
    nlevels=2^vlen
    out=fill((falses(0),UInt64(0),zeros(UInt64,0)),nlevels)
    Threads.@threads for i in 1:nlevels
        out[i]=find_matches(res_matrix,i, vlen)
    end
    return out
end

function tableCounts_dict(res_matrix::SubArray)
    tids=1:Threads.nthreads()
    resultvec=[Dict() for i in tids]
    countvec=[Dict() for i in tids]
    Threads.@threads for i in eachrow(res_matrix)
        resdict=resultvec[Threads.threadid()]
        cntdict=countvec[Threads.threadid()]
        if haskey(resdict, i)
            push!(resdict[i], i.indices[1])
            cntdict[i] += 1
        else
            resdict[i] = [i.indices[1]]
            cntdict[i] = 1
        end
    end
    finalcounts=Dict()
    finalindices=Dict()
    for i in tids
        mergewith!(append!,finalindices, resultvec[i])
        mergewith!(+, finalcounts,countvec[i])
    end
    return (finalindices, finalcounts)
end

function tableCounts_dict(res_matrix::Matrix)
    tids=1:Threads.nthreads()
    resultvec=[Dict() for i in tids]
    countvec=[Dict() for i in tids]
    Threads.@threads for i in eachrow(res_matrix)
        resdict=resultvec[Threads.threadid()]
        cntdict=countvec[Threads.threadid()]
        if haskey(resdict, i)
            push!(resdict[i], i.indices[1])
            cntdict[i] += 1
        else
            resdict[i] = [i.indices[1]]
            cntdict[i] = 1
        end
    end
    finalcounts=Dict()
    finalindices=Dict()
    for i in tids
        mergewith!(append!,finalindices, resultvec[i])
        mergewith!(+, finalcounts,countvec[i])
    end
    return (finalindices, finalcounts)
end

# TODO implement the varnames > 3 version but dict version will do for now.
function tableCounts(res_matrix,varnames)
    # if length(varnames) > 3
    return tableCounts_dict(res_matrix) 
    # else
    #return tableCounts_bf(res_matrix)
    #end
end

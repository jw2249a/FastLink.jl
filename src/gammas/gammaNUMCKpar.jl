# B-tree for difference buckets
struct IndexValue{T}
    indx::UInt32
    value::T
end

struct ValueRange{T}
    left::T
    right::T
    function ValueRange{T}(value, offset; compare! = -) where T
        new{T}(compare!(value,offset),compare!(value,-offset))
    end
end

mutable struct Bucket{T}
    range::ValueRange{T}
    obs::Vector{IndexValue{T}}
    left::Bucket{T}
    right::Bucket{T}
    lock::ReentrantLock
    function Bucket{T}(index, value, offset) where T
        n = new{T}()
        n.range = ValueRange{T}(value, offset)
        n.obs=IndexValue[IndexValue(UInt32(index), value)]
        n.lock = ReentrantLock()
        return n
    end
end

# allows for redefinition of comparison function based on need
function compare_values(ValueType::DataType,fn::Function)
    if ValueType <: Float64
        return (a::Float64,b::Float64) -> fn(a,b)
    else
        return (a::Int64,b::Int64) -> fn(a,b)
    end
end

function update_tree(cut_b,BucketImpl)
    return function update_tree!(node::Bucket,index::UInt32, value)
	if value < node.range.left
            lock(node.lock)
            if isdefined(node,:left)
                unlock(node.lock)
                update_tree!(node.left,index, value)
            else
                node.left = BucketImpl(index, value, cut_b)
                unlock(node.lock)
            end
        elseif value > node.range.right
            lock(node.lock)
            if isdefined(node,:right)
                unlock(node.lock)
                update_tree!(node.right,index, value)
            else
                node.right = BucketImpl(index, value,cut_b)
                unlock(node.lock)
            end
        else
            lock(node.lock)
            push!(node.obs, IndexValue(UInt32(index), value))
            unlock(node.lock)
        end
        
        return node
    end
end

function score_diff(results::SubArray,array_2Dindex::Function,
                partial::Bool, cut_a::Float64,cut_b::Float64,
                    match2::Vector{Bool}, match1::Vector{Bool},
                    compare!::Function)
    if partial
        (node::Bucket, value::Float64,index::UInt32) -> begin
            for i in node.obs
                if compare!(value, i.value) <= cut_a
                    results[array_2Dindex(i.indx,index),:] = match2
                elseif compare!(value, i.value) <= cut_b
                    results[array_2Dindex(i.indx,index),:] = match1
                end
            end
            return nothing
        end
    else
        (node::Bucket, value::Float64,index::UInt32) -> begin
            for i in node.obs
                if compare!(value, i.value) <= cut_a
                    results[array_2Dindex(i.indx,index),:] = match2
                end
            end
            return nothing
        end
    end
end

function score_diff(results::SubArray,array_2Dindex::Function,
                    partial::Bool, cut_a::Int,cut_b::Int,
                    match2::Vector{Bool}, match1::Vector{Bool},
                    compare!::Function)
    if partial
        (node::Bucket, value::Int,index::UInt32) -> begin
            for i in node.obs
                if abs(compare!(value, i.value)) <= cut_a
                    results[array_2Dindex(i.indx,index),:] = match2
                elseif abs(compare!(value, i.value)) <= cut_b
                    results[array_2Dindex(i.indx,index),:] = match1
                end
            end
            return nothing
        end
    else
        (node::Bucket, value::Int,index::UInt32) -> begin
            for i in node.obs
                if abs(compare!(value, i.value)) <= cut_a
                    results[array_2Dindex(i.indx,index),:] = match2
                end
            end
            return nothing
        end
    end
end

function search_upper(compare!::Function, score_diff!::Function)
    return function search_upper!(node::Bucket,index::UInt32, value, cut_b)
        diff = compare!(value, node.range.right)
        if diff > cut_b
            if isdefined(node,:right)
                search_upper!(node.right,index,value,cut_b)
            end
            # if exactly center upper will return it
        elseif diff <= -cut_b
            if isdefined(node,:left)
                search_upper!(node.left,index, value,cut_b)
            end
        else
            score_diff!(node,value,index)
            if diff > 0
                if isdefined(node,:right)
                    search_upper!(node.right,index,value,cut_b)
                end
            elseif diff < 0
                if isdefined(node,:left)
                    search_upper!(node.left,index, value,cut_b)
                end
            end
        end
        return nothing
    end
end

function search_lower(compare!::Function,score_diff!::Function)
    return function search_lower!(node::Bucket,index::UInt32, value, cut_b)
        diff = compare!(value, node.range.left)
        if diff > cut_b
            if isdefined(node,:right)
                search_lower!(node.right,index,value,cut_b)
            end
            # does not handle exactly center
        elseif diff < -cut_b
            if isdefined(node,:left)
                search_lower!(node.left,index, value,cut_b)
            end
        else
            score_diff!(node,value,index)
            if diff > 0
                if isdefined(node,:right)
                    search_lower!(node.right,index,value,cut_b)
                end
            elseif diff < 0
                if isdefined(node,:left)
                    search_lower!(node.left,index, value,cut_b)
                end
            end
        end
        return nothing
    end
end

"""
Numeric comparison of two columns

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::SubArray`: ResultMatrix object's result_matrix.
- `array_2Dindex::Function`: ResultMatrix object's array_2Dindex function
- `dims::Tuple`: ResultMatrix object's dims.
- `cut_a::Number=1.0`: Lower bound for close string distances.
- `cut_b::Number=2.0`: Lower bound for partial string distances.
"""
function gammaNUMCKpar!(results::SubArray,array_2Dindex::Function,
                        dims::Tuple; cut_a=1,cut_b=2,
                        partial::Bool=true,match2=[true,true],
                        match1=[true,false],missingval=[false,true],
                        comparison_function=-, match_method="int")


    
    ValueType=match_method=="int" ? Int64 : Float64
    compare! = compare_values(ValueType, comparison_function)    
    # define scoring function
    score_diff! = score_diff(results,
                             array_2Dindex,
                             partial,
                             cut_a,
                             cut_b,
                             match2,
                             match1,
                             compare!)

    search_upper! = search_upper(compare!,score_diff!)
    search_lower! = search_lower(compare!,score_diff!)
    
    # inner fun
    (vecA,vecB::Vector) -> begin
        # coerce to float if needed
        if typeof(vecA) <: Vector{Union{Missing, Float64}} ||
            typeof(vecB) <: Vector{Union{Missing, Float64}} ||
            match_method == "float"
            vecA=convert(Vector{Union{Missing, Float64}},vecA)
            vecB=convert(Vector{Union{Missing, Float64}},vecB)
            cut_a=convert(Float64,cut_a)
            cut_b=convert(Float64,cut_b)
            ValueType=Float64
            score_diff! = score_diff(results,
                                     array_2Dindex,
                                     partial,
                                     cut_a,
                                     cut_b,
                                     match2,
                                     match1,
                                     compare!)

            search_upper! = search_upper(compare!,score_diff!)
            search_lower! = search_lower(compare!,score_diff!)
            
        end
        
        
        
        missings_x=UInt32[]
        missings_y=UInt32[]
        missing_lock = ReentrantLock()
        # root of B-tree
        BucketImpl=Bucket{ValueType}
        update_tree! = update_tree(cut_b, BucketImpl)
        root = BucketImpl(UInt32(1), vecA[1],cut_b)

        Threads.@threads for i in UInt32(2):UInt32(length(vecA))
            value=vecA[i]
            if ismissing(value)
                lock(missing_lock)
                push!(missings_x,UInt32(i))
                unlock(missing_lock)
            else
                update_tree!(root,i, value)
            end
        end
        Threads.@threads for i in UInt32(1):UInt32(length(vecB))
            value=vecB[i]
            if ismissing(value)
                lock(missing_lock)
                push!(missings_y,i)
                unlock(missing_lock)
            else
                search_upper!(root,i, value,cut_b)
                search_lower!(root,i, value,cut_b)
            end
        end
        # set all to missing where x is missing
        if missings_x != UInt32[]
            Threads.@threads for iy in 1:dims[2]
                for ix in missings_x
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = missingval
                end
            end
        end
        # set all to missing where y is missing
        if missings_y != UInt32[]
            Threads.@threads for ix in 1:dims[1]
                for iy in missings_y
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = missingval
                end
            end
        end 

        return nothing
    end
end

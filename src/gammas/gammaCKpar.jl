# match2 = reshape(match2,(1,2))
#        match1 = reshape(match1,(1,2))
#        missingval = reshape(missingval,(1,2))
function score_strings(results::SubArray,array_2Dindex::Function,
                        distmethod::String, partial::Bool, cut_a::Float64,cut_b::Float64, w::Float64,
                        match2::Vector{Bool}, match1::Vector{Bool})

    thresh_a, thresh_b = 1 - cut_a, 1 - cut_b

    # assign distance function
    if distmethod=="jw"
        distance = JaroWinkler(p=w)
    elseif distmethod=="dl"
        distance = DamerauLevenshtein()
    elseif distmethod=="jaro"
        distance = Jaro(p=w)
    elseif distmethod=="lv"
        distance = Levenshtein()
    end
    
    if partial
        return (string1,string2,indices_x,indices_y) -> begin
            dist=round(distance(string1,string2),digits=4)
            # if matches at a threshold, go through result vector and assign new value
            if dist <= thresh_a
                for ix in indices_x,iy in indices_y
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = match2
                end
            elseif dist <= thresh_b
                for ix in indices_x,iy in indices_y
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = match1
                end
            end
            return nothing
        end
    else
        return (string1,string2,indices_x,indices_y) -> begin
            dist=round(distance(string1,string2),digits=4)
            # if matches at a threshold, go through result vector and assign new value
            if dist <= thresh_a
                for ix in indices_x,iy in indices_y
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = match2
                end
            end
            return nothing
        end
    end
end


"""
String comparison of two columns with partial match.

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::SubArray`: ResultMatrix object's result_matrix.
- `array_2Dindex::Function`: ResultMatrix object's array_2Dindex function
- `dims::Tuple`: ResultMatrix object's dims.
- `cut_a::Float=0.92`: Lower bound for close string distances.
- `cut_b::Float=0.88`: Lower bound for partial string distances.
- `distmethod::String`: String distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming).
- `w`: Winkler weight for jw string distance.
"""
function gammaCKpar!(results::SubArray,array_2Dindex::Function,dims::Tuple;
                    distmethod="jw",cut_a=0.92,cut_b=0.88,partial=true,w=0.1,
                    match2 = [true, true], match1 = [true, false],
                    missingval = [false, true])

    score_strings! = score_strings(results,array_2Dindex,distmethod,
                                    partial,cut_a,cut_b,w,match2,match1)

    return (vecA::PooledVector,vecB::PooledVector) -> begin
        # Segment unique keys from missing key
        missingvals_x = findfirst(ismissing.(vecA.pool))
        iter_x=filter(x -> x != missingvals_x, 0x00000001:UInt32(length(vecA.pool)))
        
        missingvals_y = findfirst(ismissing.(vecB.pool))
        iter_y=filter(x -> x != missingvals_y, 0x00000001:UInt32(length(vecB.pool)))
        
        # Form match matrices based on differing levels of matches
        Threads.@threads for x in iter_x
            indices_x = findall(vecA.refs .=== x)
            for y in  iter_y
                indices_y = findall(vecB.refs .=== y)
                score_strings!(vecA.pool[x],vecB.pool[y], indices_x,indices_y)
            end
        end

        # set all to missing where x is missing
        if !isnothing(missingvals_x)
            missingindices = findall(vecA.refs .== missingvals_x)
            Threads.@threads for iy in 1:dims[2]
                for ix in missingindices
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = missingval
                end
            end
        end
        # set all to missing where y is missing
        if !isnothing(missingvals_y)
            missingindices = findall(vecB.refs .== missingvals_y)
            Threads.@threads for ix in 1:dims[1]
                for iy in missingindices
                    results[array_2Dindex(UInt32(ix),UInt32(iy)),:] = missingval
                end
            end
        end
        # Return nothing
        return nothing
        
    end
end

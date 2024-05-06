function score_value2(dist::Float64,indices_x::Vector{<:Integer},indices_y::Vector{<:Integer}, cut_a::Float64, cut_b::Float64, results::DiBitMatrix)
    # if matches at a threshold, go through result vector and assign new value
    if dist >= cut_a
        for ix in indices_x, iy in indices_y
            results[ix,iy] = match2
        end
    elseif dist >= cut_b
        for ix in indices_x,iy in indices_y
            results[ix,iy] = match1
        end
    end
    return nothing
end

function score_value(dist::Float64,indices_x::Vector{<:Integer},indices_y::Vector{<:Integer}, cut_a::Float64, cut_b::Float64, results::DiBitMatrix)
    # if matches at a threshold, go through result vector and assign new value
    if dist >= cut_a
        for ix in indices_x,iy in indices_y
            results[ix,iy] = match2
        end
    end
    return nothing
end


"""
String comparison of two columns with partial match.

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::DiBitMatrix`: DiBitMatrix object's result_matrix.
- `_dims::Tuple`: DiBitMatrix object's _dims.
- `cut_a::Float=0.92`: Lower bound for close string distances.
- `cut_b::Float=0.88`: Lower bound for partial string distances.
- `distmethod::String`: String distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming).
- `w`: Winkler weight for jw string distance.
"""
function gammaCKpar!(vecA::PooledVector,vecB::PooledVector,
                     results::DiBitMatrix;
                     distmethod="jw",cut_a=0.92,cut_b=0.88,partial=true,w=0.1)

    if @isdefined(_dims) == false
        _dims = (length(vecA), length(vecB))
    end
    
    # assign distance function
    if distmethod=="jw"
        distance = JaroWinkler(p=w)
    elseif distmethod=="dl"
        distance = DamerauLevenshtein()
    elseif distmethod=="jaro"
        distance = Jaro(p=w)
    elseif distmethod=="lv"
        distance = Levenshtein()
    elseif distmethod=="ro"
        distance = RatcliffObershelp()
    elseif distmethod=="osa"
        distance = OptimalStringAlignment()
    elseif distmethod=="hamming"
        distance = Hamming()
    end


    if partial
        score_value! = score_value2
    else
        score_value! = score_value
    end
    
    # Segment unique keys from missing key
    missingvals_x = findfirst(ismissing.(vecA.pool))
    iter_x=filter(x -> x != missingvals_x, UInt32(1):UInt32(length(vecA.pool)))
    
    missingvals_y = findfirst(ismissing.(vecB.pool))
    iter_y=filter(x -> x != missingvals_y, UInt32(1):UInt32(length(vecB.pool)))
    
    # Form match matrices based on differing levels of matches
    Threads.@threads for x in iter_x
        indices_x = findall(vecA.refs .=== x)
        for y in  iter_y
            indices_y = findall(vecB.refs .=== y)
            dist=round(compare(vecA.pool[x],vecB.pool[y], distance),digits=4) #this always normalizes dist 0 to 1
            score_value!(dist, indices_x,indices_y, cut_a,cut_b, results)
        end
    end

    # set all to missing where x is missing
    if !isnothing(missingvals_x)
        missingindices = findall(vecA.refs .== missingvals_x)
        Threads.@threads for iy in 1:_dims[2]
            for ix in missingindices
                results[ix,iy] = missingval
            end
        end
    end
    
    # set all to missing where y is missing
    if !isnothing(missingvals_y)
        missingindices = findall(vecB.refs .== missingvals_y)
        Threads.@threads for ix in 1:_dims[1]
            for iy in missingindices
                results[ix,iy] = missingval
            end
        end
    end
    # Return nothing
    return nothing
end

# term frequency adjusted version
function gammaCKpar!(vecA::PooledVector,vecB::PooledVector,
                     results::DiBitMatrix,
                     tf_table_x::SubArray{Float16},
                     tf_table_y::SubArray{Float16};
                     distmethod="jw",cut_a=0.92,cut_b=0.88,partial=true,w=0.1,
                     tf_minimum_u_value=0.001)

    if @isdefined(_dims) == false
        _dims = (length(vecA), length(vecB))
    end
    
    
    # assign distance function
    if distmethod=="jw"
        distance = JaroWinkler(p=w)
    elseif distmethod=="dl"
        distance = DamerauLevenshtein()
    elseif distmethod=="jaro"
        distance = Jaro(p=w)
    elseif distmethod=="lv"
        distance = Levenshtein()
    elseif distmethod=="ro"
        distance = RatcliffObershelp()
    elseif distmethod=="osa"
        distance = OptimalStringAlignment()
    elseif distmethod=="hamming"
        distance = Hamming()
    end

    if partial
        score_value! = score_value2
    else
        score_value! = score_value
    end
    
    # Segment unique keys from missing key
    missingvals_x = findfirst(ismissing.(vecA.pool))
    iter_x=filter(x -> x != missingvals_x, UInt32(1):UInt32(length(vecA.pool)))
    
    missingvals_y = findfirst(ismissing.(vecB.pool))
    iter_y=filter(x -> x != missingvals_y, UInt32(1):UInt32(length(vecB.pool)))
    
    # Form match matrices based on differing levels of matches
    Threads.@threads for x in iter_x
        # all values in x that match unique value
        indices_x = findall(vecA.refs .=== x)
        
        # term frequency adjustment for x
        tf_val_x = length(indices_x)/_dims[1]
        for tf_i in indices_x
            tf_table_x[tf_i] = max(tf_val_x, tf_minimum_u_value)
        end
        
        for y in iter_y
            # all values in y that match unique value
            indices_y = findall(vecB.refs .=== y)
            
            # term frequency adjustment for y
            tf_val_y = length(indices_y)/_dims[2]
            for tf_i in indices_y
                tf_table_y[tf_i] = max(tf_val_y, tf_minimum_u_value)
            end

            # string comparison
            dist=round(compare(vecA.pool[x],vecB.pool[y], distance),digits=4) #this always normalizes dist 0 to 1
            score_value!(dist, indices_x,indices_y, cut_a,cut_b, results)
        end
    end

    # set all to missing where x is missing
    if !isnothing(missingvals_x)
        missingindices = findall(vecA.refs .== missingvals_x)

        # term frequency adjustment for x
        tf_val_x = length(missingindices)/_dims[1]
        for tf_i in missingindices
            tf_table_x[tf_i] = max(tf_val_y, tf_minimum_u_value)
        end
        
        Threads.@threads for iy in 1:_dims[2]
            for ix in missingindices
                results[ix,iy] = missingval
            end
        end
    end
    
    # set all to missing where y is missing
    if !isnothing(missingvals_y)
        missingindices = findall(vecB.refs .== missingvals_y)
        # term frequency adjustment for y
        tf_val_y = length(missingindices)/_dims[2]
        for tf_i in missingindices
            tf_table_y[tf_i] = max(tf_val_y, tf_minimum_u_value)
        end

        Threads.@threads for ix in 1:_dims[1]
            for iy in missingindices
                results[ix,iy] = missingval
            end
        end
    end
    # Return nothing
    return nothing
end


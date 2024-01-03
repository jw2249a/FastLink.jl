function gammaCKpar!(vecA::PooledVector,vecB::PooledVector, results::SubArray, array_2Dindex::Function, dims::Tuple;
                     cut_a=0.92, cut_b=0.88, distmethod="jw", w=0.1)
    match2 = [true, true]
    match1 = [true, false]
    missingval = [false, true]
    # assign distance function
    if distmethod=="jw"
        distance = JaroWinkler(p=w)
    end
    
    # Threshold vals are 1 - cuts
    thresh_a, thresh_b = 1 - cut_a, 1 - cut_b

    # Segment unique keys from missing key
    missingvals_x=findfirst(ismissing.(vecA.pool))
    iter_x=filter(x -> x != missingvals_x, 0x00000001:UInt32(length(vecA.pool)))
    
    missingvals_y=findfirst(ismissing.(vecB.pool))
    iter_y=filter(x -> x != missingvals_y, 0x00000001:UInt32(length(vecB.pool)))
    
    # Form match matrices based on differing levels of matches
    Threads.@threads for x in iter_x
        indices_x = findall(vecA.refs .=== x)
        for y in  iter_y
            indices_y = findall(vecB.refs .=== y)
            dist=round(distance(vecA.pool[x],vecB.pool[y]),digits=4)
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

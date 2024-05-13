"""
Exact comparison of variables. 

# Arguments
- `vecA::PooledVector`: Target column of dfB for string comparison.
- `vecB::PooledVector`: Target column of dfB for string comparison.
- `results::SubArray`: ResultMatrix object's result_matrix.
- `array_2Dindex::Function`: ResultMatrix object's array_2Dindex function
- `_dims::Tuple`: ResultMatrix object's _dims.
"""
function gammaKpar!(vecA::PooledVector,vecB::PooledVector,results::DiBitMatrix)
    if @isdefined(_dims) == false
        _dims = (length(vecA), length(vecB))
    end
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
            # if matches at a threshold, go through result vector and assign new value
            if vecA.pool[x] == vecB.pool[y]
                for ix in indices_x,iy in indices_y
                    results[ix,iy] = match2
                end
            end
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
function gammaKpar!(vecA::PooledVector,vecB::PooledVector,results::DiBitMatrix,
                    tf_table_x::SubArray{Float16},
                    tf_table_y::SubArray{Float16};
                    tf_minimum_u_value=0.001)

    if @isdefined(_dims) == false
        _dims = (length(vecA), length(vecB))
    end
    # Segment unique keys from missing key
    missingvals_x = findfirst(ismissing.(vecA.pool))
    iter_x=filter(x -> x != missingvals_x, 0x00000001:UInt32(length(vecA.pool)))
    
    missingvals_y = findfirst(ismissing.(vecB.pool))
    iter_y=filter(x -> x != missingvals_y, 0x00000001:UInt32(length(vecB.pool)))
    
    # Form match matrices based on differing levels of matches
    Threads.@threads for x in iter_x
        indices_x = findall(vecA.refs .=== x)
         # term frequency adjustment for x
        tf_val_x = length(indices_x)/_dims[1]
        for tf_i in indices_x
            tf_table_x[tf_i] = max(tf_val_x, tf_minimum_u_value)
        end
        for y in  iter_y
            indices_y = findall(vecB.refs .=== y)
             # term frequency adjustment for y
            tf_val_y = length(indices_y)/_dims[2]
            for tf_i in indices_y
                tf_table_y[tf_i] = max(tf_val_y, tf_minimum_u_value)
            end
            # if matches at a threshold, go through result vector and assign new value
            if vecA.pool[x] == vecB.pool[y]
                for ix in indices_x,iy in indices_y
                    results[ix,iy] = match2
                end
            end
        end
    end

    # set all to missing where x is missing
    if !isnothing(missingvals_x)
        missingindices = findall(vecA.refs .== missingvals_x)
        # term frequency adjustment for x
        tf_val_x = length(missingindices)/_dims[1]
        for tf_i in missingindices
            tf_table_x[tf_i] = max(tf_val_x, tf_minimum_u_value)
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

# vector version of gammakpar for highly entropic data (aka high ordinality compared to obs)
function gammaKpar!(vecA::Vector,vecB::Vector,results::DiBitMatrix)
    if @isdefined(_dims) == false
        _dims = (length(vecA), length(vecB))
    end
    # Segment unique keys from missing key
    missingvals_x = findfirst(ismissing.(vecA))
    iter_x=filter(x -> x != missingvals_x, 0x00000001:UInt32(length(vecA)))
    
    missingvals_y = findfirst(ismissing.(vecB))
    iter_y=filter(x -> x != missingvals_y, 0x00000001:UInt32(length(vecB)))
    
    # Form match matrices based on differing levels of matches
    Threads.@threads for (ix, x) in collect(enumerate(vecA))
        indices_x = findall(vecA .=== x)
        if ismissing(x)
            for iy in collect(1:dims[2])
                results[ix,iy] = missingval
            end
        else
            Threads.@threads for (iy, y) in collect(enumerate(vecB))
                if ismissing(y)
                    results[ix,iy] = missingval
                elseif x == y
                    results[ix,iy] = match2
                end
                # if matches at a threshold, go through result vector and assign new value
            end
        end
    end
    # Return nothing
    return nothing
end


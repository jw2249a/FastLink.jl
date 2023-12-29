# needed to create ranges for the result matrix
function create_ranges(vec)
    counter=1
    weight=0
    ranges = []
    for i in vec
        push!(ranges, counter:(i+weight))
        counter += i
        weight += i
    end
    return ranges
end

# prespecify matrix that is able to hold results of comparisons
struct ResultMatrix
    comparison_bits::Vector{Int}
    dims::Tuple{Int, Int}
    ranges::Vector{UnitRange{Int}}
    result_matrix::Matrix{Bool}
    array_2Dindex::Function
    
    function ResultMatrix(comparison_bits::Vector{Int}, dims::Tuple{Int, Int})
        comparison_bits=Int.(ceil.(log2.(comparison_bits.+2)))
        new(comparison_bits, dims,
            create_ranges(comparison_bits),
            falses((prod(dims),sum(comparison_bits))),
            function array_2Dindex(row::Int, col::Int; nrows::Int = dims[1])
                return row + (col - 1) * nrows
            end)
    end
end



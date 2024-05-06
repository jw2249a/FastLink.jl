module matchpatterns
import Base: getindex, setindex!
import DataFrames: DataFrame

using ..DiBitMat

export ComparisonIndex, LocalPatterns, MatchPatterns

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

# extending base with get and set based on comparison index
getindex(df::DataFrame, x::ComparisonIndex) = df[x.row,x.col]
setindex!(df::DataFrame, value, x::ComparisonIndex) = df[x.row,x.col] = value

function getindex(vm::DiBitMatrix, x::ComparisonIndex)
    linear_index = (x.col - 1) * vm.nrows + x.row
    return vm.data[linear_index]
end
function setindex!(vm::DiBitMatrix, value::UInt8, x::ComparisonIndex)
        linear_index = (x.col - 1) * vm.nrows + x.row
        vm.data[linear_index] = value
end


end # module MatchPatterns



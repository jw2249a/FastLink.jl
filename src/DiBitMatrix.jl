module DiBitMat
import Base: getindex, setindex!, view
import DataStructures: DiBitVector
export DiBitMatrix

"""
Extending DiBitVectors from DataStructures.jl to include matrices.
"""
struct DiBitMatrix
    data::DiBitVector
    nrows::Integer
    ncols::Integer
end

# base definition of the DiBitMatrix
function DiBitMatrix(nrows::Integer, ncols::Integer)
    data = DiBitVector(nrows * ncols, 0)  # Or choose an appropriate type
    return DiBitMatrix(data, nrows, ncols)
end

# getting items by index
function getindex(vm::DiBitMatrix, i::Int, j::Int)
    linear_index = (j - 1) * vm.nrows + i
    return vm.data[linear_index]
end

function getindex(vm::DiBitMatrix, ::Colon, j::Int)
    column = zeros(UInt8, vm.nrows)
    for i in 1:vm.nrows
        linear_index = (j - 1) * vm.nrows + i
        column[i] = vm.data[linear_index]
    end
    return column
end

function getindex(vm::DiBitMatrix, i::Int, ::Colon)
    row = zeros(UInt8, vm.ncols)
    for j in 1:vm.ncols
        linear_index = (j - 1) * vm.nrows + i
        row[j] = vm.data[linear_index]
    end
    return row
end

# setting items by index
function setindex!(vm::DiBitMatrix, value::UInt8, i::T, j::T) where {T<:Integer}
    linear_index = (j - 1) * vm.nrows + i
    vm.data[linear_index] = value
end

# extending view to handle DiBitMatrix columns
function getIndices(vm::DiBitMatrix,::Colon,j::Int)
    return (j - 1) * vm.nrows + 1, (j - 1) * vm.nrows + vm.nrows
end

function getIndices(vm::DiBitMatrix, i::Int, ::Colon)
    row = zeros(Integer, vm.ncols)
    for j in 1:vm.ncols
        row[j] = (j - 1) * vm.nrows + i
    end
    return row
end

function view(vm::DiBitMatrix,::Colon, j::Int)
    start,finish=getIndices(vm,:,j)
    return view(vm.data, start:finish)
end

function view(vm::DiBitMatrix, i::Int,::Colon)
    vals=getIndices(vm, i,:)
    return view(vm.data, vals)
end


end

module FastLink
using DataFrames
import PooledArrays: PooledVector
import Distributions: Dirichlet,rand

# match constants
const match1::UInt8 = UInt8(1)
const match2::UInt8 = UInt8(2)
const missingval::UInt8 = UInt8(3)

include("DiBitMatrix.jl")
using .DiBitMat
include("gammas/Gammas.jl")
using .Gammas

include("matchPatterns.jl")
include("emlink.jl")
include("getMatches.jl")
include("fastlink/fastlink.jl")

export(fastLink)


end # module FastLink

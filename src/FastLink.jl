module FastLink
using DataFrames
import PooledArrays: PooledVector
import Distributions: Dirichlet,rand


include("resultMatrix.jl")
include("gammas/Gammas.jl")
using .Gammas

include("tableCounts.jl")
include("emlink.jl")
include("getMatches.jl")
include("fastlink/fastlink.jl")

export(tableCounts)
export(fastLink)


end # module FastLink

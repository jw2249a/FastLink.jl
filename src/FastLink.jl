module FastLink
using DataFrames
using PooledArrays
using Distributions
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein


include("resultMatrix.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")
include("gammaKpar.jl")
include("tableCounts.jl")
include("emlink.jl")
include("getMatches.jl")
include("utils/prettyprinting.jl")
include("fastlink/fastlink.jl")

export(tableCounts)
export(fastLink)


end # module FastLink

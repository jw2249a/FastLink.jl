module FastLink
using DataFrames
using PooledArrays
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein


include("resultMatrix.jl")
include("gammaCKpar.jl")
include("fastlink/fastlink.jl")

end # module FastLink

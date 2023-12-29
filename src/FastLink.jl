module FastLink
using DataFrames
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein
import PooledArrays: PooledVector

include("resultMatrix.jl")
include("gammaCKpar.jl")
include("fastlink/fastlink.jl")

end # module FastLink

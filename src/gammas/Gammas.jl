module Gammas
import PooledArrays: PooledVector, PooledArray
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein, compare

using ..DiBitMat
import ..nonmatch, ..match1, ..match2, ..missingval


include("gammaKpar.jl")
include("gammaNUMCKpar.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")
export gammaCKpar!, gammaKpar!, gammaCKfuzzy!, gammaNUMCKpar!

end # Gammas

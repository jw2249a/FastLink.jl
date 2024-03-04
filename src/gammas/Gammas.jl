module Gammas
import PooledArrays: PooledVector
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein, compare
using FastLink.DiBitMat
import FastLink: match1, match2, missingval, nonmatch

include("gammaKpar.jl")
include("gammaNUMCKpar.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")
export gammaCKpar!, gammaKpar!, gammaCKfuzzy!, gammaNUMCKpar!

end # Gammas

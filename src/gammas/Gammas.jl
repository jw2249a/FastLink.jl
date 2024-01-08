module Gammas
import PooledArrays: PooledVector
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein

include("gammaKpar.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")

export gammaCKpar!, gammaKpar!, gammaCKfuzzy!

end # Gammas

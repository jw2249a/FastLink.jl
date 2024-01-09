module Gammas
import PooledArrays: PooledVector
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein

include("gammaKpar.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")
include("gammaNUMCKpar.jl")
export gammaCKpar!, gammaKpar!, gammaCKfuzzy!, gammaNUMCKpar!

end # Gammas

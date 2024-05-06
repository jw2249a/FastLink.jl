module FastLink
using DataFrames
import PooledArrays: PooledVector

# match constants
const nonmatch::UInt8 = UInt8(0)
const match1::UInt8 = UInt8(1)
const match2::UInt8 = UInt8(2)
const missingval::UInt8 = UInt8(3)

const STRING_DISTANCE_METHODS = Dict("jw" => "jw",
                                     "jarowinkler" => "jw",
                                     "jaro winkler" => "jw",
                                     "jaro-winkler" => "jw",
                                     "jaro" => "jaro",
                                     "dl" => "dl",
                                     "dameraulevenshtein" => "dl",
                                     "damerau levenshtein" => "dl",
                                     "damerau-levenshtein" => "dl",
                                     "lv" => "lv",
                                     "levenshtein" => "lv",
                                     "hamming" => "hamming",
                                     "ro" => "ro",
                                     "ratcliffobershelp" => "ro",
                                     "ratcliff obershelp" => "ro",
                                     "ratcliff-obershelp" => "ro",
                                     "osa" => "osa",
                                     "optimal string alignment" => "osa",
                                     "optimalstringalignment" => "osa"
                                     )

include("settings/settings.jl")
include("DiBitMatrix.jl")
include("matchPatterns.jl")
include("gammas/Gammas.jl")
include("term_frequency_adjustment.jl")
include("emlink.jl")
include("patterns.jl")

using .settings
using .DiBitMat
using .matchpatterns
using .Gammas
using .emlink
using .tf
using .patterns

include("fastlink/fastlink.jl")

export gammaCKpar!, gammaKpar!, gammaCKfuzzy!, gammaNUMCKpar!, DiBitMatrix, namedtuple, fetch_parameters, retrieve, parse_configuration, remove_keys, emlinkMARmov,  STRING_DISTANCE_METHODS, match1, match2, missingval, nonmatch, indices_to_uids, process_comparisons, fastLink

#export(fastLink)


end # module FastLink

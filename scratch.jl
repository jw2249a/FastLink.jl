using DataFrames
using BenchmarkTools
using CSV

#varnames=["FIRST_NAME"]
#varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME"]
varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]

N1=10
N2=10
include("utils/prettyprinting.jl")


using DataFrames
using PooledArrays
import StringDistances: Jaro, JaroWinkler, Levenshtein, DamerauLevenshtein


include("resultMatrix.jl")
include("gammaCKpar.jl")
include("gammaCKfuzzy.jl")
include("tableCounts.jl")
include("utils/prettyprinting.jl")
include("fastlink/fastlink.jl")

a_fil="../../../rstudio/test_merge/data/test_a.csv"
b_fil="../../../rstudio/test_merge/data/test_b.csv"
dfA=CSV.read(a_fil, DataFrame,
             limit=N1,
             ignoreemptyrows=true,
             ntasks=1,
             pool=true,
             missingstring=["", "NA", "NaN", "NULL", "Null"])

dfB=CSV.read(b_fil, DataFrame,
             limit=N2,
             ignoreemptyrows=true,
             ntasks=1,
             pool=true,
             missingstring=["", "NA", "NaN", "NULL", "Null"])

allowmissing!(dfA)
allowmissing!(dfB)



function reverse_array_2Dindex(index, nrows)::Tuple{UInt32, UInt32}
    # Adjust index to 0-based for calculation
    zero_based_index = index - 1

    # Calculate row and column
    col = UInt32(div(zero_based_index, nrows)) + 1
    row = UInt32(mod(zero_based_index, nrows)) + 1

    return (row, col)
end

function get_matches(vecA, vecB, results; match_pattern=[1, 1])
    targs=findall(x->x == match_pattern,eachrow(results))
    return [(vecA[i],vecB[x]) for (i,x) in reverse_array_2Dindex.(targs,length(vecA))]
end

vecA=dfB.FIRST_NAME
vecB=dfB.FIRST_NAME
vars=fastLinkVars(dfA, dfB, varnames)
comparison_levels=[2 for i in varnames]
res=ResultMatrix(comparison_levels, (length(vecA), length(vecB)))

#res=ResultMatrix(comparison_levels, (nrow(dfA), nrow(dfB)))
#results=view(res.result_matrix,:,res.ranges[1])

gammaCKfuzzy!(vecA,
              vecB,
              view(res.result_matrix,:,res.ranges[1]),
              res.array_2Dindex,
              res.dims)
[(i,ii) for (i,ii) in get_matches(vecA,vecB,res.result_matrix[:,1:2],match_pattern=[1,1])]

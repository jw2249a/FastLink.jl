using Pkg
Pkg.develop(path=".")
Pkg.precompile()
using DataFrames
using BenchmarkTools
using CSV
using FastLink
using PooledArrays
import Pkg.Artifacts: @artifact_str


a_fil = @artifact_str "dfA"
b_fil = @artifact_str "dfB"

varnames=["firstname","middlename", "lastname","housenum"]
match_method=["string", "string","string", "float"]
cut_a=[0.92,0.92,0.92,1]
cut_p=[0.88,0.88,0.88,2]



dfA=CSV.read("$(a_fil)/dfA.csv", DataFrame,
             ntasks=1,
             pool=true,
             missingstring=["", "NA"])
dfB=CSV.read("$(b_fil)/dfB.csv", DataFrame,
             ntasks=1,
             pool=true,
             missingstring=["", "NA"])
dfA.id = hash.(eachrow(dfA))
dfB.id2 = hash.(eachrow(dfB))


for var in varnames[1:3]
    dfA[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfA[:,var]))
    dfB[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfB[:,var]))
end

results=fastLink(dfA,dfB,varnames,("id","id2"),
                 match_method=match_method,
                 term_freq_adjustment=[true],
                 cut_a=cut_a,cut_p=cut_p,
                 threshold_match = 0.85)





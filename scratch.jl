using FastLink
using DataFrames
using CSV
using PooledArrays
import Pkg.Artifacts: @artifact_str
#Pkg.add(url="https://github.com/jw2249a/FastLink.jl")
using StatsBase
using JSON

a_fil = @artifact_str "dfA"
b_fil = @artifact_str "dfB"


dfA=CSV.read("$(a_fil)/dfA.csv", DataFrame,
             ntasks=1,
             pool=true,
             missingstring=["", "NA"])
dfB=CSV.read("$(b_fil)/dfB.csv", DataFrame,
             ntasks=1,
             pool=true,
             missingstring=["", "NA"])


config = JSON.parsefile("test_parameters.json")

dfA.id = hash.(eachrow(dfA))
dfB.id2 = hash.(eachrow(dfB))


varnames=["firstname","middlename", "lastname","housenum"]

for var in varnames
    if eltype(dfA[:,var]) <: AbstractString
        dfA[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfA[:,var]))
        dfB[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfB[:,var]))
    else
        dfA[!,var] = Vector(dfA[!,var])
        dfB[!,var] = Vector(dfB[!,var])
    end
end

result=fastLink(dfA, dfB, config)


using DataFrames
using FastLink
using BenchmarkTools
using CSV

N1=10_000
N2=100_000
fil=""
vf_fil=""


vf=CSV.read(vf_fil, DataFrame,
            limit=N1,
            ignoreemptyrows=true,
            ntasks=1,
            pool=true,
            missingstring=["", "NA", "NaN", "NULL", "Null"])

tv=CSV.read(fil, DataFrame,
            limit=N2,
            ignoreemptyrows=true,
            ntasks=1,
            pool=true,
            missingstring=["", "NA", "NaN", "NULL", "Null"])

varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME"]


@time FastLink.fastLink(tv,vf,varnames)






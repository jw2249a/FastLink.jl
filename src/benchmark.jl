using DataFrames
using FastLink
using BenchmarkTools
using CSV

N1=10_000
N2=10_000
fil="../../rstudio/test_merge/data/tv_export_CA.csv"
vf_fil="../../rstudio/test_merge/data/vf_export_CA.csv"


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



comparison_levels=[2,1,2,1,2]
res=ResultMatrix(comparison_levels, (nrow(tv), nrow(vf)))

@btime FastLink.fastLink(tv,vf,varnames) evals=1

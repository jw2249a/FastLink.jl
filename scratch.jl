using Pkg
Pkg.develop(path=".")
using DataFrames
using BenchmarkTools
using CSV
using FastLink
using PooledArrays

numeric=true
# files for performance
test=false
if test
    a_fil="../dfA.csv"
    b_fil="../dfB.csv"
    if numeric
        varnames=["housenum"]
        match_method=["float"]
        cut_a=[1]
        cut_p=[2]
    else
        varnames=["firstname","middlename", "lastname","streetname"]
        match_method=["string", "string","string", "string"]
        cut_a=[0.92,0.92,0.92,0.92]
        cut_p=[0.88,0.88,0.88,0.88]
    end
else
    a_fil="../../rstudio/test_merge/data/test_a.csv"
    b_fil="../../rstudio/test_merge/data/test_b.csv"

    if numeric
        varnames=["ZIP", "DOB_YEAR", "ZIP4"]
        match_method=["float", "float", "float"]
        cut_a=[1,1,1]
        cut_p=[2,2,2]
    else
        varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME"]
        cut_a=[0.92,0.92,0.92,0.92]
        cut_p=[0.88,0.88,0.88,0.88]
        #varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]
    end
end


#[100,200,500,1_000,2_000,4_000, 5_000, 10_000,20_000, 40_000, 50_000,100_000,1_000_000]
N1=1_000
N2=100_000


if test
    dfA=CSV.read(a_fil, DataFrame,
                 ntasks=1,
                 pool=true,
                 missingstring=["", "NA"])
    dfB=CSV.read(b_fil, DataFrame,
                 ntasks=1,
                 pool=true,
                 missingstring=["", "NA"])
else
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
end


if !test && numeric
    for var in varnames
        dfA[!,var]=passmissing(x-> try return parse(Float64,x) catch e return 0.0 end).(dfA[:,var])
        dfB[!,var]=passmissing(x-> try return parse(Float64,x) catch e return 0.0 end).(dfB[:,var])
    end
end

if test && !numeric
    for var in varnames
        dfA[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfA[:,var]))
        dfB[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfB[:,var]))
    end
end



@btime results=fastLink(dfA,dfB,varnames,match_method=match_method,cut_a=cut_a,cut_p=cut_p)()



using DataFrames
using BenchmarkTools
using CSV
using FastLink
using PooledArrays

# files for performance
test=true
if test
    a_fil="../../rstudio/test_merge/data/dfA.csv"
    b_fil="../../rstudio/test_merge/data/dfB.csv"
    varnames=["firstname","middlename", "lastname","streetname"]
    else
    a_fil="../../rstudio/test_merge/data/test_a.csv"
    b_fil="../../rstudio/test_merge/data/test_b.csv"
    varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME"]
    #varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]
end


#[100,200,500,1_000,2_000,4_000, 5_000, 10_000,20_000, 40_000, 50_000,100_000,1_000_000]
N1=1_000
N2=1_000


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


if test
    for var in varnames
        dfA[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfA[:,var]))
        dfB[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfB[:,var]))
    end
end


results=fastLink(dfA,dfB,varnames,false)


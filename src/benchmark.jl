using DataFrames
using FastLink
using BenchmarkTools
using CSV

include("utils/prettyprinting.jl")

a_fil="../../rstudio/test_merge/data/test_a.csv"
b_fil="../../rstudio/test_merge/data/test_b.csv"

#varnames=["FIRST_NAME"]
#varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME"]
varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]

N1=1_000
N2_N=[100,200,500,1_000,2_000,4_000, 5_000, 10_000,20_000, 40_000, 50_000,100_000,1_000_000]
println("## $(length(varnames)) vars")
for N2 in N2_N

    println(center_in_line("( $(pretty_print_number(N1)) x $(pretty_print_number(N2)) )"))
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

    
    GC.gc()
    @btime FastLink.fastLink($dfA,$dfB,$varnames)
    println("")
end




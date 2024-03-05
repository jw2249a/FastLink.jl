using DataFrames
using FastLink
using BenchmarkTools
using CSV
import Pkg.Artifacts: @artifact_str
using Profile

include("utils/prettyprinting.jl")

a_fil="../../rstudio/test_merge/data/test_a.csv"
b_fil="../../rstudio/test_merge/data/test_b.csv"

#varnames=["FIRST_NAME"]
varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME"]
match_type=["fuzzy","fuzzy","fuzzy","fuzzy"]
#varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]
#[100,200,500,1_000,2_000,4_000, 5_000, 10_000,20_000, 40_000, 50_000,100_000,1_000_000]
idvars=("TS_ID","TV_ID")
N2=20_000
N1_N=[10_000,50_000,100_000,500_000,750_000,1_000_000]
println("## $(length(varnames)) vars")
for N1 in N1_N

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

    @btime fastLink($dfA,$dfB,$varnames,$idvars,match_method=$match_type)

end

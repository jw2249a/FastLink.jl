using Base: postoutput
using Pkg
Pkg.develop(path="..")
using FastLink
using DataFrames
using BenchmarkTools
using CSV
import Pkg.Artifacts: @artifact_str
using Profile

outputfile = "../benchmark.csv"

include("utils/prettyprinting.jl")

a_fil="../../../rstudio/test_merge/data/test_a.csv"
b_fil="../../../rstudio/test_merge/data/test_b.csv"


config = Dict("link_type"=>"link_only",
     "idvar"=> ["TV_ID", "TS_ID"],
     "comparisons"=> Dict("name" => "total",
                         "threshold_match" => 0.88,
                         "variables" => [
                             Dict("varname" => "FIRST_NAME",
                                  "method" => "jarowinkler",
                                  "tf_adjust" => true),
                             Dict("varname" => "MIDDLE_NAME",
                                  "method" => "exact",
                                  "tf_adjust" => true),
                             Dict("varname" => "STREET_NAME",
                                  "method" => "jarowinkler",
                                  "tf_adjust" => true) ]))


config_tf = Dict("link_type"=>"link_only",
     "idvar"=> ["TV_ID", "TS_ID"],
     "comparisons"=> Dict("name" => "total",
                         "threshold_match" => 0.88,
                         "variables" => [
                             Dict("varname" => "FIRST_NAME",
                                  "method" => "jarowinkler"),
                             Dict("varname" => "MIDDLE_NAME",
                                  "method" => "exact"),
                             Dict("varname" => "STREET_NAME",
                                  "method" => "jarowinkler")]))

open(outputfile, "w" do file
         write(file,"\"N1\",\"N2\",\"u_FIRST_NAME\",\"u_MIDDLE_NAME\"")


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
    
    res = @benchmark fastLink($dfA,$dfB,$config)

end

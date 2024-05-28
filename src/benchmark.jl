using Pkg
Pkg.develop(path="..")
using FastLink
using DataFrames
using BenchmarkTools
using CSV
import Pkg.Artifacts: @artifact_str


outputfile = "../benchmark.csv"
outputfile_tf = "../benchmark_tf.csv"
include("utils/prettyprinting.jl")

a_fil="../../../rstudio/test_merge/data/test_a.csv"
b_fil="../../../rstudio/test_merge/data/test_b.csv"

sims=10
N1=10_000
N2_N=[100,200,500,1000,2000,5000,10_000,20_000,50_000,100_000,200_000,500_000,750_000,1_000_000]


config = Dict("link_type"=>"link_only",
     "idvar"=> ["TV_ID", "TS_ID"],
     "comparisons"=> Dict("name" => "total",
                         "threshold_match" => 0.88,
                         "variables" => [
                             Dict("varname" => "FIRST_NAME",
                                  "method" => "jarowinkler"),
                             Dict("varname" => "FIRST_NAME1",
                                  "method" => "jarowinkler"),
                             Dict("varname" => "MIDDLE_NAME",
                                  "method" => "exact"),
                             Dict("varname" => "MIDDLE_NAME1",
                                  "method" => "exact"),
                             Dict("varname" => "STREET_NAME",
                                  "method" => "jarowinkler"),
                             Dict("varname" => "STREET_NAME1",
                                  "method" => "jarowinkler")]))

config_tf = Dict("link_type"=>"link_only",
     "idvar"=> ["TV_ID", "TS_ID"],
     "comparisons"=> Dict("name" => "total",
                         "threshold_match" => 0.88,
                         "variables" => [
                             Dict("varname" => "FIRST_NAME",
                                  "method" => "jarowinkler",
                                  "tf_adjust" => true),
                             Dict("varname" => "FIRST_NAME1",
                                  "method" => "jarowinkler",
                                  "tf_adjust" => true),
                             Dict("varname" => "MIDDLE_NAME",
                                  "method" => "exact",
                                  "tf_adjust" => true),
                              Dict("varname" => "MIDDLE_NAME1",
                                  "method" => "exact",
                                  "tf_adjust" => true),
                             Dict("varname" => "STREET_NAME",
                                  "method" => "jarowinkler",
                                  "tf_adjust" => true),
                          Dict("varname" => "STREET_NAME1",
                               "method" => "jarowinkler",
                               "tf_adjust" => true) ]))

# creating output files for the benchmarks
vars=retrieve(config, "varname")
filheadernames=append!(["allocs", "time_elapsed", "N_dfA", "N_dfB", "sim_num"],
                       "uniqueobsa_" .* vars,
                       "uniqueobsb_" .* vars,
                       "time_" .* vars)

filheader=*((v != last(filheadernames) ? "\"$v\"," : "\"$v\"\n" for v in filheadernames)...)
open(outputfile, "w") do file
    write(file, filheader)
end
open(outputfile_tf, "w") do file
    write(file, filheader)
end


for sim_num in 1:sims
    for N2 in N2_N

        @info center_in_line("(sim $sim_num, $(pretty_print_number(N1)) x $(pretty_print_number(N2)) )")
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

        dfA.STREET_NAME1 = Vector(dfA.STREET_NAME)
        dfB.STREET_NAME1 = Vector(dfB.STREET_NAME)
        dfA.FIRST_NAME1 = Vector(dfA.FIRST_NAME)
        dfB.FIRST_NAME1 = Vector(dfB.FIRST_NAME)
        dfA.MIDDLE_NAME1 = Vector(dfA.MIDDLE_NAME)
        dfB.MIDDLE_NAME1 = Vector(dfB.MIDDLE_NAME)
        @info "running no tf"
        let res_allocs = @allocated fastLink(dfA,dfB,config) 
            
            stime = time()
            res = fastLink(dfA,dfB,config, true)
            telapsed = time() - stime

            res=DataFrame(append!(["allocs"=>res_allocs, "time_elapsed"=>telapsed,"N1"=>N1,"N2"=>N2, "sim_num"=>sim_num],
                                  ["uniqueobsa_$v"=>length(unique(dfA[:,v])) for v in vars],
                                  ["uniqueobsb_$v"=>length(unique(dfB[:,v])) for v in vars],
                                  ["time_$(vars[i])"=>res["benchtimes"][i] for i in 1:length(vars)]))
            CSV.write(outputfile, res, writeheader = false, append = true)
        end
        @info "running tf"
        let res_allocs = @allocated fastLink(dfA,dfB,config) 
            
            stime = time()
            res = fastLink(dfA,dfB,config_tf, true)
            telapsed = time() - stime

            res=DataFrame(append!(["allocs"=>res_allocs,"time_elapsed"=>telapsed,"N1"=>N1,"N2"=>N2, "sim_num"=>sim_num],
                                  ["uniqueobsa_$v"=>length(unique(dfA[:,v])) for v in vars],
                                  ["uniqueobsb_$v"=>length(unique(dfB[:,v])) for v in vars],
                                  ["time_$(vars[i])"=>res["benchtimes"][i] for i in 1:length(vars)]))
            CSV.write(outputfile_tf, res, writeheader = false, append = true)
        end
        
        
        
    end
end

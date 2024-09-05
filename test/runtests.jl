using FastLink
using Test
import DataFrames: DataFrame, nrow, passmissing, PooledArray
import CSV
import Pkg.Artifacts: @artifact_str


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

dfA.ida = hash.(eachrow(dfA))
dfB.idb = hash.(eachrow(dfB))

varnames=["firstname","middlename", "lastname","housenum"]

for var in varnames[1:3]
    dfA[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfA[:,var]))
    dfB[!,var] = PooledArray(passmissing(x->uppercase(x)).(dfB[:,var]))
end


dfA.housenum = Vector(dfA.housenum)
dfB.housenum = Vector(dfB.housenum)
config = Dict("idvar" => ["ida", "idb"],
              "link_type" => "link_only",
              "comparisons" => Dict("name" => "total",
                                    "variables" => [
                                        Dict("varname" => "firstname",
                                             "partial" => true,
                                             "method" => "jarowinkler",
                                             "cut_a" => 0.92,
                                             "cut_b" => 0.88,
                                             "w" => 0.1),
                                        Dict("varname" => "middlename",
                                             "partial" => true,
                                             "method" => "jarowinkler",
                                             "cut_a" => 0.92,
                                             "cut_b" => 0.88,
                                             "w" => 0.1),
                                        Dict("varname" => "lastname",
                                             "partial" => true,
                                             "method" => "jarowinkler",
                                             "cut_a" => 0.92,
                                             "cut_b" => 0.88,
                                             "w" => 0.1),
                                        Dict("varname" => "housenum",
                                             "partial" => true,
                                             "method" => "numeric",
                                             "cut_a" => 1,
                                             "cut_b" => 2)
                                    ]))

@testset "Testing FastLink Basic Run" begin
    @info "Executing fastLink()"
    results=fastLink(dfA,dfB,config)
    @info "Correct # of Matches"
    matches = getMatches(results)
    p_w = results["resultsEM"]["patterns_w"]
    inds = p_w.zeta_j .>= results["resultsEM"]["threshold_match"]
    @test sum(p_w.counts[inds]) == 50
    @info "Correct grouping of matches"
    @info p_w.counts[inds] == nrow.(matches)
    @info "Number of patterns == 26"
    @test results["resultsEM"]["number_of_unique_patterns"] == 26
    @info "Number of counts == (N₁×N₂)"
    @test sum(p_w.counts) == nrow(dfA) * nrow(dfB)
    @info "Ρ(𝑢) >=.999"
    @test results["resultsEM"]["p_u"] >= .999
    @info "Ρ(𝑚) <= .0005"
    @test results["resultsEM"]["p_m"] <= .0005
    true
end

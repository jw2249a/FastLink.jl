
using FastLink
using Test
import DataFrames: DataFrame
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

varnames = ["firstname","middlename", "lastname","housenum"]
cut_a = [0.92,0.92,0.92,1]
cut_p = [0.88,0.88,0.88,2]
match_method = ["string","string","string","float"]
partials = [true,true,true,true]
fuzzy = [false,false,false,false]
stringdist_method = ["jw","jw","jw","jw"]
upper_case = [false,false,false,false]
jw_weight = [0.1,0.1,0.1,0.0]




@testset "Testing FastLink Basic Run" begin
    @info "Executing fastLink()"
    results=fastLink(dfA,dfB,varnames,
                     match_method=match_method,
                     partials=partials,
                     fuzzy=fuzzy,
                     upper_case=upper_case,
                     stringdist_method=stringdist_method,
                     cut_a=cut_a,
                     cut_p=cut_p,
                     jw_weight=jw_weight)()
    @test true

    @info "Correct # of Matches"
    @test length(results[2]) == 50
    @info "Number of patterns == 26"
    @test length(results[1].patterns_b) == 26
    @info "Number of counts == (N₁×N₂) "
    @test sum(results[1].patterns_w.counts) == nrow(dfA) * nrow(dfB)
    @info "Ρ(𝑢) >=.999"
    @test results[1].p_u >= .999
    @info "Ρ(𝑚) <= .0005"
    @test results[1].p_m <= .0005


    
    
    true
end


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

@testset "FastLink.jl" begin
    varnames = ["firstname","middlename", "lastname","housenum"]
    cut_a = [0.92,0.92,0.92,1]
    cut_p = [0.88,0.88,0.88,2]
    match_method = ["string","string","string","float"]
    partials = [true,true,false,true]
    fuzzy = [true,false,true,false]
    stringdist_method = ["jw","jw","jw",""]
    upper_case = [false,false,false,false]
    jw_weight = [0.1,0.1,0.1,0.0]
    
    fastLink(dfA,dfB,varnames,
             match_method=match_method,
             partials=partials,
             fuzzy=fuzzy,
             upper_case=upper_case,
             stringdist_method=stringdist_method,
             cut_a=cut_a,
             cut_p=cut_p,
             jw_weight=jw_weight)()

    println("completed")
    return true
end

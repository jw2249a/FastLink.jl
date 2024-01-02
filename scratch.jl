using DataFrames
using BenchmarkTools
using CSV
using PooledArrays
using Distributions
using StringDistances
include("src/resultMatrix.jl")
include("src/tableCounts.jl")
include("src/gammaCKfuzzy.jl")
include("src/gammaCKpar.jl")
include("src/fastlink/fastlink.jl")
include("src/utils/prettyprinting.jl")

# files for performance
test=true
if test
    a_fil="../../rstudio/test_merge/data/dfA.csv"
    b_fil="../../rstudio/test_merge/data/dfB.csv"
    varnames=["firstname", "lastname"]
    else
    a_fil="../../rstudio/test_merge/data/test_a.csv"
    b_fil="../../rstudio/test_merge/data/test_b.csv"
    varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME"]
    #varnames=["FIRST_NAME", "MIDDLE_NAME", "LAST_NAME", "STREET_NAME", "STATE"]
end



#varnames=["FIRST_NAME"]
comparison_levels=[2 for i in varnames]

#[100,200,500,1_000,2_000,4_000, 5_000, 10_000,20_000, 40_000, 50_000,100_000,1_000_000]
N2=1_000
N1=1_000

nrow(dfB)

if test
    dfA=CSV.read(a_fil, DataFrame,
                 ntasks=1,
                 pool=true)

    dfB=CSV.read(b_fil, DataFrame,
                 ntasks=1,
                 pool=true)

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

function get_match_list(N)
    return (results[N][1],Int(results[N][2]))
end
function reverse_array_2Dindex(index::Int,var,dfA::DataFrame,dfB::DataFrame)
    nrows= nrow(dfA)
    # Adjust index to 0-based for calculation
    zero_based_index = index - 1

    # Calculate row and column
    row = Int(mod(zero_based_index, nrows)) + 1
    col = Int(div(zero_based_index, nrows)) + 1


    return ([row,col],(dfA[row,var], dfB[col,var]))
end


res.result_matrix[:,1:4]



results2=findall(x-> x== [true,false,false,false],eachrow(res.result_matrix[:,1:4])) .|> Tuple .|> first .|> x->reverse_array_2Dindex(x,["firstname","lastname"],dfA,dfB)



filter(y->y[1][1]==43,results2)

results2[524880]

[get_match_list(i) for i in 1:16] 


missingval = [false, true]
match2 = [true,true]
match1 = [true,false]

obs_a = nrow(dfA)
obs_b = nrow(dfB)
res=ResultMatrix(comparison_levels, (obs_a, obs_b))

tol = Float64(1e-05)

delta = Float64(1)
count = 1


patterns = deepcopy(results)
nfeatures=length(varnames)


𝛾ⱼₖ=collect(keys(patterns[2]))
nⱼ = collect(values(patterns[2]))
N = length(nⱼ)
pₘ = 0.1



# if λ priors are not declared
ψ = 1
μ = 1

# if π priors are null
## for address
⍺₀_address = 1
⍺₁_address = 1
address_field = falses(nfeatures)
## for lambda 
⍺₀_gender = 1
⍺₁_gender = 1
genderaddress_field = falses(nfeatures)

pᵤ = 1 - pₘ

# 𝛾ₖ(𝑖,𝑗) = m ~idp Discrete(πₖₘ)
# Initializing matrices for conditional probabilities
prob_𝛾𝑘𝑚 = fill(Float64[], nfeatures)
prob_𝛾𝑘𝑢 = fill(Float64[], nfeatures)
vals_𝛾𝑗𝑘 = fill(Vector{Union{Missing,Int}}(),nfeatures)
uvals_𝛾𝑗𝑘 = fill(Int[],nfeatures)

for c in 1:nfeatures
    col=res.ranges[c]
    vals_𝛾𝑗𝑘[c] = [i[col] == missingval ? missing : sum(i[col]) for i in 𝛾ⱼₖ]
    uvals_𝛾𝑗𝑘[c] = sort(unique([i for i in vals_𝛾𝑗𝑘[c] if !ismissing(i)]))
    c_m = collect(1:50:(length(uvals_𝛾𝑗𝑘[c])*50))
    prob_𝛾𝑘𝑚[c] = sort(rand(Dirichlet(c_m),1)[:],rev=true)
    prob_𝛾𝑘𝑢[c] = sort(rand(Dirichlet(c_m),1)[:],rev=false)
end

prob_𝛾𝑘𝑗𝑚 = missings(Union{Missing,Float64}, (nfeatures,N))
prob_𝛾𝑘𝑗𝑢 = missings(Union{Missing,Float64}, (nfeatures,N))


prob_𝛾𝑗𝑢 = missings(Union{Missing,Float64}, N)


delta >= tol

for i in 1:nfeatures
    prob_𝛾𝑘𝑗𝑚[i,:] = [ismissing(j) ? j : prob_𝛾𝑘𝑚[i][findfirst(uvals_𝛾𝑗𝑘[i] .== j)] for j in vals_𝛾𝑗𝑘[i]]
    prob_𝛾𝑘𝑗𝑢[i,:] = [ismissing(j) ? j : prob_𝛾𝑘𝑢[i][findfirst(uvals_𝛾𝑗𝑘[i] .== j)] for j in vals_𝛾𝑗𝑘[i]]
end

prob_𝛾𝑗𝑚 = sum.(skipmissing.(eachcol(log.(prob_𝛾𝑘𝑗𝑚))))
prob_𝛾𝑗𝑢 = sum.(skipmissing.(eachcol(log.(prob_𝛾𝑘𝑗𝑢))))

log_prod_𝛾𝑗𝑚 = prob_𝛾𝑗𝑚 .+ log(pₘ) 
log_prod_𝛾𝑗𝑢 = prob_𝛾𝑗𝑢 .+ log(pᵤ)

logxpy = function(lx,ly)
    return maximum.(eachrow([lx ly])) .+ log1p.(exp.(-abs.(lx .- ly)))
end

probability_vector = function(x::Vector{Float64})
    return x./sum(x)
end


𝜁ⱼ = exp.(log_prod_𝛾𝑗𝑚 - logxpy(log_prod_𝛾𝑗𝑚,log_prod_𝛾𝑗𝑢))
num_prod = exp.(log.(nⱼ) + log.(𝜁ⱼ))

pₘ = exp(log(sum(num_prod) + μ - 1) - log(ψ - μ + sum(nⱼ)))
pᵤ = 1-pₘ


# for i in 1:nfeatures
prob_𝛾𝑘𝑚[i]

vals_𝛾𝑗𝑘[i]

num_prod

ismissing.(vals_𝛾𝑗𝑘[i]) .== false 

vals_𝛾𝑗𝑘[i].==uvals_𝛾𝑗𝑘[i][1]

prob_𝛾𝑘𝑚[i,:] = 
[sum(num_prod[findall(ismissing.(vals_𝛾𝑗𝑘[i]) .== false .&& vals_𝛾𝑗𝑘[i].==uvals_𝛾𝑗𝑘[i][j])]) for j in 1:length(uvals_𝛾𝑗𝑘[i])]

probability_vector([let sub1= sub=findall(ismissing.(vals_𝛾𝑗𝑘[i]) .== false .&& vals_𝛾𝑗𝑘[i].==uvals_𝛾𝑗𝑘[i][j]);
     sum(nⱼ[sub] - num_prod[sub]) end for j in 1:length(uvals_𝛾𝑗𝑘[i])])

skipmissing(vals_𝛾𝑗𝑘[i])

[ismissing(j) ? j : prob_𝛾𝑘𝑚[i][findfirst(uvals_𝛾𝑗𝑘[i] .== j)] for j in vals_𝛾𝑗𝑘[i]]
[ismissing(j) ? j : prob_𝛾𝑘𝑢[i][findfirst(uvals_𝛾𝑗𝑘[i] .== j)] for j in vals_𝛾𝑗𝑘[i]]

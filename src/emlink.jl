function recursive_flatten(x)
    result = []
    for item in x
        if item isa AbstractArray
            append!(result, recursive_flatten(item))
        else
            push!(result, item)
        end
    end
    return result
end
function logxpy(lx,ly)
    return maximum.(eachrow([lx ly])) .+ log1p.(exp.(-abs.(lx .- ly)))
end

function probability_vector(x::Vector{Float64})
    return x./sum(x)
end 

function skipmissing_equality(x, y)
    return  replace(x .== y, missing=>false)
end

"""
Expectation maximization function. 
"""
function emlinkMARmov(patterns::Dict, obs_a::Int,obs_b::Int,varnames::Vector{String}, ranges::Vector{UnitRange{Int64}}; p_m=0.1,iter_max=5000,tol=Float64(1e-05),missingval = [false,true])
    # Initialize count and delta for while loop and break point
    delta = Float64(1)
    count = 1
    
    # Info for EM algorithm
    p_u = 1 - p_m
    nfeatures=length(varnames)
    gamma_jk=collect(keys(patterns))
    n_j = collect(values(patterns))
    N = length(n_j)
   
    # TODO: add "if statement" λ priors are declared
    psi = 1
    mu = 1
    
    ###########################################
    # # TODO: add "if statement" for π priors #
    # ## for address                          #
    # ⍺₀_address = 1                          #
    # ⍺₁_address = 1                          #
    # address_field = falses(nfeatures)       #
    # ## for lambda                           #
    # ⍺₀_gender = 1                           #
    # ⍺₁_gender = 1                           #
    # genderaddress_field = falses(nfeatures) #
    ###########################################
    # initialize variables that need value to be returned
    zeta_j=0.0
    num_prod = zeros(Float64,0)
    p_old=zeros(Float64, 0)
    p_gamma_jm = zeros(Union{Missing,Float64},0)
    p_gamma_ju = zeros(Union{Missing,Float64},0)
    # Initializing matrices for conditional probabilities
    p_gamma_km = fill(Float64[], nfeatures)
    p_gamma_ku = fill(Float64[], nfeatures)
    vals_gamma_jk = fill(Vector{Union{Missing,Int}}(),nfeatures)
    uvals_gamma_jk = fill(Int[],nfeatures)
    p_gamma_kjm = missings(Union{Missing,Float64}, (nfeatures,N))
    p_gamma_kju = missings(Union{Missing,Float64}, (nfeatures,N))
    for c in 1:nfeatures
        col=ranges[c]
        vals_gamma_jk[c] = [i[col] == missingval ? missing : sum(i[col]) for i in gamma_jk]
        uvals_gamma_jk[c] = sort(unique([i for i in vals_gamma_jk[c] if !ismissing(i)]))
        c_m = collect(1:50:(length(uvals_gamma_jk[c])*50))
        p_gamma_km[c] = sort(rand(Dirichlet(c_m),1)[:],rev=false)
        p_gamma_ku[c] = sort(rand(Dirichlet(c_m),1)[:],rev=true)
    end
    
    while abs(delta) >= tol
        p_old=recursive_flatten([p_m,p_u,p_gamma_km,p_gamma_ku])
        for i in 1:nfeatures
            p_gamma_kjm[i,:] = [ismissing(j) ? missing : p_gamma_km[i][findfirst(uvals_gamma_jk[i] .== j)] for j in vals_gamma_jk[i]]
            p_gamma_kju[i,:] = [ismissing(j) ? missing : p_gamma_ku[i][findfirst(uvals_gamma_jk[i] .== j)] for j in vals_gamma_jk[i]]
        end
        p_gamma_jm = sum.(skipmissing.(eachcol(log.(p_gamma_kjm))))
        p_gamma_ju = sum.(skipmissing.(eachcol(log.(p_gamma_kju))))
        log_prod_gamma_jm = p_gamma_jm .+ log(p_m)
        log_prod_gamma_ju = p_gamma_ju .+ log(p_u)
        zeta_j = exp.(log_prod_gamma_jm - logxpy(log_prod_gamma_jm,log_prod_gamma_ju))
        num_prod = exp.(log.(n_j) + log.(zeta_j))
        p_m = exp(log(sum(num_prod) + mu - 1) - log(psi - mu + sum(n_j)))
        p_u = 1-p_m

        for i in 1:nfeatures
            p_gamma_km[i] = 
                sort(probability_vector([sum(num_prod[findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]))])
                      for j in 1:length(uvals_gamma_jk[i])]),rev=false)
            p_gamma_ku[i] =
                sort(probability_vector([let sub1 = sub=findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]));
                                             sum(n_j[sub] - num_prod[sub])
                                         end
                                         for j in 1:length(uvals_gamma_jk[i])]),rev=true)
        end
        delta = maximum(abs.(recursive_flatten([p_m,p_u,p_gamma_km,p_gamma_ku]) - p_old))
        
        count +=1
        if count > iter_max
            @warn "The EM algorithm has run for the specified number of iterations but has not converged yet."
            break
        end
    end
    weights = p_gamma_jm - p_gamma_ju
    p_gamma_jm=probability_vector(exp.(p_gamma_jm))
    p_gamma_ju=probability_vector(exp.(p_gamma_ju))
    data_gamma_ = DataFrame(vals_gamma_jk, :auto)
    rename!(data_gamma_,["gamma_$i" for i in 1:ncol(data_gamma_)])
    data_w=hcat(data_gamma_,DataFrame((counts=n_j,weights,p_gamma_jm,p_gamma_ju)))
    
    return (zeta_j = zeta_j, p_m = p_m, p_u = p_u,
            pgamma_km = p_gamma_km, pgamma_ku = p_gamma_ku,
            pgamma_jm = p_gamma_jm, pgamma_ju = p_gamma_ju,
            patterns_w = data_w,
            patterns_b = gamma_jk,
            iter_converge = count,
            obs_a = obs_a, obs_b = obs_b,
            varnames = varnames)
end

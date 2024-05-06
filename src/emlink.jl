module emlink
using DataFrames
import Distributions: Dirichlet,rand
using ..matchpatterns
import ..nonmatch, ..match1, ..match2, ..missingval
export emlinkMARmov

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

function probability_vector(x::Vector{BigFloat})
    return x./sum(x)
end 

function skipmissing_equality(x, y)
    return  replace(x .== y, missing=>false)
end



"""
Expectation maximization function. 
"""
function emlinkMARmov(patterns::MatchPatterns, dims::Tuple{Int,Int},varnames::Vector{String};
                      p_m=0.1,iter_max=5000,tol=1e-05, prior_lambda=0.0, w_lambda=0.0,
                      prior_pi=0.0, w_pi=0.0, address_field=Vector{Bool}(), threshold_match=0.85, u_b=1e10,
                      digit_precision=16)
    obs_a, obs_b = dims

    # precision for comparisons
    old=precision(BigFloat)
    setprecision(BigFloat,digit_precision; base=10)
    atexit(x->setprecision(BigFloat,old))
    # Initialize count and delta for while loop and break point
    delta = BigFloat(1)
    count = 1
    p_m=BigFloat(p_m)
    # Info for EM algorithm
    p_u = 1 - p_m
    nfeatures=length(varnames)
    gamma_jk=patterns.patterns
    n_j = length.(patterns.indices)
    N = length(n_j)
    psi = 0
    mu = 0
    alpha0_address = 1.0
    alpha1_address = 1.0
    # if λ priors are declared
    if prior_lambda == 0
        psi = 1
        mu = 1        
    else
        if w_lambda == 0
            @error "If declaring a lambda prior, you need to declare weights via w_lambda."
        elseif w_lambda > 1 || w_lambda < 0
            @error "w_lambda must be between 0 and 1."
        elseif w_lambda == 1
            w_lambda = 1 - 1e-05
        end
        c_lambda = w_lambda/(1-w_lambda)
        # hyperparameters for lambda
        mu = prior_lambda * c_lambda * obs_a * obs_b + 1
        psi = (1 - prior_lambda) * mu / prior_lambda
    end
    
    # if pi prior is declared
    if prior_pi == 0
        alpha0_address = 1
        alpha1_address = 1
        address_field = falses(nfeatures)
    else
        if prior_lambda == 0
            @error "If declaring a prior on pi, you need to declare lambda prior."
        elseif w_pi == 0
            @error "If providing a prior for pi, please specify the weight using w_pi"
        elseif w_pi < 0 | w_pi > 1
            @error "w_pi must be between 0 and 1."
        elseif w_pi == 1
            w_pi = 1 - 1e-05
        end

        c_pi = w_pi / (1 - w_pi)
        exp_match = prior_lambda * obs_a * obs_b

        # Optimal hyperparameters for pi
        alpha0_address = c_pi * prior_pi * exp_match + 1
        alpha1_address = alpha0_address * (1 - prior_pi) / prior_pi
    end

    # initialize variables that need value to be returned
    zeta_j=BigFloat(0)
    num_prod = zeros(BigFloat,0)
    p_old=zeros(BigFloat, 0)
    p_gamma_jm = zeros(Union{Missing,BigFloat},0)
    p_gamma_ju = zeros(Union{Missing,BigFloat},0)
    # Initializing matrices for conditional probabilities
    p_gamma_km = fill(BigFloat[], nfeatures)
    p_gamma_ku = fill(BigFloat[], nfeatures)
    vals_gamma_jk = fill(Vector{Union{Missing,Int}}(),nfeatures)
    uvals_gamma_jk = fill(Int[],nfeatures)
    p_gamma_kjm = missings(Union{Missing,BigFloat}, (nfeatures,N))
    p_gamma_kju = missings(Union{Missing,BigFloat}, (nfeatures,N))
    for c in 1:nfeatures
        vals_gamma_jk[c] = [i[c] == missingval ? missing : sum(i[c]) for i in gamma_jk]
        uvals_gamma_jk[c] = sort(unique([i for i in vals_gamma_jk[c] if !ismissing(i)]))
        c_m = BigFloat.(collect(1:50:(length(uvals_gamma_jk[c])*50)))
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
            km_prob=sort([sum(num_prod[findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]))])
                          for j in 1:length(uvals_gamma_jk[i])],rev=false)
            if address_field[i]
                km_prob += append!([alpha0_address], [alpha1_address for i in 1:(length(uvals_gamma_jk[i])-1)])
            end
            p_gamma_km[i] = 
                probability_vector(km_prob)
            p_gamma_ku[i] =
                sort(probability_vector([let sub1 = sub=findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]));
                                             sum(n_j[sub] - num_prod[sub])
                                         end
                                         for j in 1:length(uvals_gamma_jk[i])]),rev=true)
        end
        delta = maximum(BigFloat.(abs.(recursive_flatten([p_m,p_u,p_gamma_km,p_gamma_ku]) - p_old)))
        
        count += 1
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

    ismatch = zeta_j .>= threshold_match .&& data_w.weights .<= u_b

    data_w.zeta_j = zeta_j
    
    
    return Dict("p_m" => p_m, "p_u" => p_u,
                "number_of_comparisons" => *(dims...),
                "number_of_unique_patterns" => nrow(data_w),
                "iter_converge" => count,
                "threshold_match" => threshold_match,
                "varnames" => varnames, 
                "patterns_w" => data_w,
                "pgamma_km" => p_gamma_km, "pgamma_ku" => p_gamma_ku,
                "pgamma_jm" => p_gamma_jm, "pgamma_ju" => p_gamma_ju)
end



function emlinkMARmov(gamma_jk::Vector{Vector{UInt8}},n_j::Vector{Int64}, dims::Tuple{Int,Int},varnames::Vector{String};
                      p_m=0.1,iter_max=5000,tol=1e-05, prior_lambda=0.0, w_lambda=0.0,
                      prior_pi=0.0,w_pi=0.0, address_field=Vector{Bool}(), threshold_match=0.85, u_b=1e10,
                      digit_precision=16)
    obs_a, obs_b = dims

    # precision for comparisons
    old=precision(BigFloat)
    setprecision(BigFloat,digit_precision; base=10)
    atexit(x->setprecision(BigFloat,old))
    # Initialize count and delta for while loop and break point
    delta = BigFloat(1)
    count = 1
    p_m=BigFloat(p_m)
    # Info for EM algorithm
    p_u = 1 - p_m
    nfeatures=length(varnames)
    N = length(n_j)
    psi = 0
    mu = 0
    alpha0_address = 1.0
    alpha1_address = 1.0
    # if λ priors are declared
    if prior_lambda == 0
        psi = 1
        mu = 1        
    else
        if w_lambda == 0
            @error "If declaring a lambda prior, you need to declare weights via w_lambda."
        elseif w_lambda > 1 || w_lambda < 0
            @error "w_lambda must be between 0 and 1."
        elseif w_lambda == 1
            w_lambda = 1 - 1e-05
        end
        c_lambda = w_lambda/(1-w_lambda)
        # hyperparameters for lambda
        mu = prior_lambda * c_lambda * obs_a * obs_b + 1
        psi = (1 - prior_lambda) * mu / prior_lambda
    end
    
    # if pi prior is declared
    if prior_pi == 0
        alpha0_address = 1
        alpha1_address = 1
        address_field = falses(nfeatures)
    else
        if prior_lambda == 0
            @error "If declaring a prior on pi, you need to declare lambda prior."
        elseif w_pi == 0
            @error "If providing a prior for pi, please specify the weight using w_pi"
        elseif w_pi < 0 | w_pi > 1
            @error "w_pi must be between 0 and 1."
        elseif w_pi == 1
            w_pi = 1 - 1e-05
        end

        c_pi = w_pi / (1 - w_pi)
        exp_match = prior_lambda * obs_a * obs_b

        # Optimal hyperparameters for pi
        alpha0_address = c_pi * prior_pi * exp_match + 1
        alpha1_address = alpha0_address * (1 - prior_pi) / prior_pi
    end

    # initialize variables that need value to be returned
    zeta_j=BigFloat(0)
    num_prod = zeros(BigFloat,0)
    p_old=zeros(BigFloat, 0)
    p_gamma_jm = zeros(Union{Missing,BigFloat},0)
    p_gamma_ju = zeros(Union{Missing,BigFloat},0)
    # Initializing matrices for conditional probabilities
    p_gamma_km = fill(BigFloat[], nfeatures)
    p_gamma_ku = fill(BigFloat[], nfeatures)
    vals_gamma_jk = fill(Vector{Union{Missing,Int}}(),nfeatures)
    uvals_gamma_jk = fill(Int[],nfeatures)
    p_gamma_kjm = missings(Union{Missing,BigFloat}, (nfeatures,N))
    p_gamma_kju = missings(Union{Missing,BigFloat}, (nfeatures,N))
    for c in 1:nfeatures
        vals_gamma_jk[c] = [i[c] == missingval ? missing : sum(i[c]) for i in gamma_jk]
        uvals_gamma_jk[c] = sort(unique([i for i in vals_gamma_jk[c] if !ismissing(i)]))
        c_m = BigFloat.(collect(1:50:(length(uvals_gamma_jk[c])*50)))
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
            km_prob=sort([sum(num_prod[findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]))])
                          for j in 1:length(uvals_gamma_jk[i])],rev=false)
            if address_field[i]
                km_prob += append!([alpha0_address], [alpha1_address for i in 1:(length(uvals_gamma_jk[i])-1)])
            end
            p_gamma_km[i] = 
                probability_vector(km_prob)
            p_gamma_ku[i] =
                sort(probability_vector([let sub1 = sub=findall(skipmissing_equality(vals_gamma_jk[i], uvals_gamma_jk[i][j]));
                                             sum(n_j[sub] - num_prod[sub])
                                         end
                                         for j in 1:length(uvals_gamma_jk[i])]),rev=true)
        end
        delta = maximum(BigFloat.(abs.(recursive_flatten([p_m,p_u,p_gamma_km,p_gamma_ku]) - p_old)))
        
        count += 1
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
    
    ismatch = zeta_j .>= threshold_match .&& data_w.weights .<= u_b
    
    data_w.zeta_j = zeta_j
    
    return Dict("p_m" => p_m, "p_u" => p_u,
                "number_of_comparisons" => *(dims...),
                "number_of_unique_patterns" => nrow(data_w),
                "iter_converge" => count,
                "threshold_match" => threshold_match,
                "varnames" => varnames, 
                "patterns_w" => data_w,
                "pgamma_km" => p_gamma_km, "pgamma_ku" => p_gamma_ku,
                "pgamma_jm" => p_gamma_jm, "pgamma_ju" => p_gamma_ju)
end

end

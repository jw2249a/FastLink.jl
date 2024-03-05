function check_input_lengths(x, var_length::Int, varname::String)
    if typeof(x) == "string"
        x=[x]
    end
    if length(x) == var_length
        return x
    elseif length(x) == 1
        return [x[1] for i in 1:var_length]
    else
        @error "Number of inputs for $varname is > 1 and < $var_length (the number of vars declared)."
    end
end

# Ensure that vartypes match each other. if string then both must be strings
# TODO: implement tighter parameter checks and coercion where needed. 
function check_var_types(x::DataFrame, y::DataFrame, varnames::Vector{String},match_method::Vector{String}, partials::Vector{Bool})
    xtypes=eltype.(eachcol(select(x, varnames)))
    ytypes=eltype.(eachcol(select(y, varnames)))
    # if partials is empty then all partial
    if partials == []
        partials = [true for i in xtypes]
    end
    comparison_levels = Int[]
    # defaults if no match methods declared
    if match_method == []
        for (ix,iy,iv,partial) in zip(xtypes,ytypes, varnames, partials)
            comparison_level = typeof(ix) <: Union || typeof(iy) <: Union || partial ? 2 : 1
            if (ix <: Union{Missing,AbstractString}) && (iy <: Union{Missing,AbstractString})
                match_type="string"
            elseif (ix <: Union{Missing,Number}) && (iy <: Union{Missing,Number})
                if ix <: Union{Missing, Float64} || iy <: Union{Missing,Float64}
                    match_type="float"
                else
                    match_type="numeric"
                end
            elseif (ix <: Union{Missing,Bool}) && (iy <: Union{Missing,Bool})
                match_type="bool"
            else
                @error "*(VAR $iv)*: dfA type $ix does not match dfB type $ix or type not known for matching"
            end
            push!(match_method, match_type)
            push!(comparison_levels, comparison_level)
        end
    else 
        for (ix,iy,iv,im,partial) in zip(xtypes,ytypes, varnames,match_method, partials)
            comparison_level = typeof(ix) <: Union || typeof(iy) <: Union || partial ? 2 : 1
            push!(comparison_levels, comparison_level)
        end
    end
    return match_method,comparison_levels
end

function remove_no_matched_var_indices(resultsEM)
    for i in eachindex(resultsEM.patterns_b)
        if (match1 in resultsEM.patterns_b[i]) ⊽ (match2 in resultsEM.patterns_b[i])
            resultsEM.indices[i] = Vector{ComparisonIndex}()
        end
    end
end


"""
Probabilistic record matching using FastLink data-matching algorithm.
Algorithm taken from:
- Enamorado, Ted, Benjamin Fifield, and Kosuke Imai. 2017. fastLink: Fast Probabilistic Record Linkage with Missing Data. Version 0.6.

# Arguments
- `dfA::DataFrame`: The first dataset to be matched.
- `dfB::DataFrame`: The second dataset to be matched.
- `varnames::Vector{String}`: Variable names to be used in the matching.
- `fuzzy::Bool`: Whether to match using a fuzzy string distance for speed (default false).
- `string_case::String`: "upper" or "lower" (only if fuzzy)
- `stringdist_method::String`: String distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming).
- `cut_a::Float`: Upper bound for string distance cutoff.
- `cut_p::Float`: Lower bound for string distance (if varnames in partial).
- `jw_weight`: Winkler weight for jw string distance.
- `tol_em`: Convergence tolerance for the EM Algorithm. (default 1e-04)
- `threshold_match`: Lower bound for the posterior probability that will act as a cutoff for matches.
- `dedupe_matches`: Whether to dedupe the matches within the dataset.

# Returns
- `MatchedData::DataFrame`: The resulting DataFrame after matching.

# Examples
```julia
matched_data = fastLink(dfA, dfB, ["firstname", "lastname", "city"])
"""
function fastLink(dfA::DataFrame, dfB::DataFrame,
                  varnames::Vector{String},
                  idvar::Tuple{String,String};
                  match_method=String[],
                  partials=[true],
                  upper_case=[true],
                  stringdist_method = ["jw"],
                  cut_a = [0.92], cut_p = [0.88],
                  jw_weight = [0.1],
                  address_field = [false],
                  tol_em = 1e-05,
                  prior_lambda = 0.0,
                  w_lambda = 0.0,
                  prior_pi = 0.0,
                  w_pi = 0.0,
                  threshold_match = 0.85,
                  dedupe_matches = true,
                  verbose = false)
    # dims
    numvars=length(varnames)
    obs_a=nrow(dfA)
    obs_b=nrow(dfB)
    dims = (obs_a,obs_b)
    
    @info "Checking settings for $numvars declared variables."
    partials = check_input_lengths(partials, numvars, "partials")
    upper_case = check_input_lengths(upper_case, numvars, "upper_case")
    jw_weight = check_input_lengths(jw_weight, numvars, "jw_weight")
    cut_a = check_input_lengths(cut_a, numvars, "cut_a")
    cut_p = check_input_lengths(cut_p, numvars, "cut_p")
    address_field = check_input_lengths(address_field, numvars, "address_field")
    stringdist_method = check_input_lengths(stringdist_method, numvars, "stringdist_method")
    
    vartypes, comparison_levels = check_var_types(dfA,dfB,varnames,match_method,partials)
    
    res = [DiBitMatrix(obs_a,obs_b) for _ in varnames]
    
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)

    
    # iterate through variables and execute function over them
    for i in eachindex(varnames)
        @info "Now matching var $(varnames[i]) using $(match_method[i])"
        if match_method[i] == "fuzzy"
            gammaCKfuzzy!(dfA[!,varnames[i]],
                          dfB[!,varnames[i]],
                          res[i],
                          dims,
                          cut_a=cut_a[i], 
                          cut_b=cut_p[i],
                          upper=upper_case[i],
                          w=jw_weight[i],
                          partial=partials[i])
        elseif match_method[i] == "string"
            gammaCKpar!(dfA[!,varnames[i]],
                        dfB[!,varnames[i]],
                        res[i],
                        dims,
                        distmethod=stringdist_method[i],
                        cut_a=cut_a[i], 
                        cut_b=cut_p[i],
                        w=jw_weight[i],
                        partial=partials[i])
        elseif match_method[i] == "exact" || match_method[i] == "bool"
            gammaKpar!(dfA[!,varnames[i]],
                       dfB[!,varnames[i]],
                       res[i],
                       dims)
        elseif match_method == "numeric" || match_method=="float" || match_method == "int"
            gammaNUMCKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           cut_a=cut_a[i],
                           cut_b=cut_p[i],
                           partial=partials[i])
        end
    end
    @info "Getting table counts"
    counts = get_match_patterns(res)

    @info "Running expectation maximization function"
    resultsEM = emlinkMARmov(counts, obs_a,obs_b,
                             varnames,tol=tol_em,
                             prior_lambda=prior_lambda, w_lambda=w_lambda,
                             prior_pi=prior_pi, w_pi=w_pi,
                             address_field=address_field)
    
    # testing removing uncessessary indices (where no obs exist)
    #remove_no_matched_var_indices(resultsEM)
    # adding uids 
    resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),))
    
    @info "Retrieving matches"
    getMatches(resultsEM,threshold_match=threshold_match)
    
    return (resultsEM) 
 
end

# version of fast link that i can pass to via named tuples
function fastLink(dfA::DataFrame, dfB::DataFrame;
                  varnames=String[],
                  match_method=String[],
                  idvar=String[],
                  partials=[true],
                  upper_case=[true],
                  stringdist_method = ["jw"],
                  cut_a = [0.92], cut_p = [0.88],
                  jw_weight = [0.1],
                  address_field = [false],
                  tol_em = 1e-05,
                  prior_lambda = 0.0,
                  w_lambda = 0.0,
                  prior_pi = 0.0,
                  w_pi = 0.0,
                  threshold_match = 0.85,
                  dedupe_matches = true,
                  verbose = false)

    # idvar to Tuple
    idvar = (idvar[1],idvar[2])
    # dims
    numvars=length(varnames)
    obs_a=nrow(dfA)
    obs_b=nrow(dfB)
    dims = (obs_a,obs_b)
    
    @info "Checking settings for $numvars declared variables."
    partials = check_input_lengths(partials, numvars, "partials")
    upper_case = check_input_lengths(upper_case, numvars, "upper_case")
    jw_weight = check_input_lengths(jw_weight, numvars, "jw_weight")
    cut_a = check_input_lengths(cut_a, numvars, "cut_a")
    cut_p = check_input_lengths(cut_p, numvars, "cut_p")
    stringdist_method = check_input_lengths(stringdist_method, numvars, "stringdist_method")
    
    vartypes, comparison_levels = check_var_types(dfA,dfB,varnames,match_method,partials)
    
    res = [DiBitMatrix(obs_a,obs_b) for _ in varnames]
    
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)

    
    # iterate through variables and execute function over them
    for i in eachindex(varnames)
        @info "Now matching var $(varnames[i]) using $(match_method[i])"
        if match_method[i] == "fuzzy"
            gammaCKfuzzy!(dfA[!,varnames[i]],
                          dfB[!,varnames[i]],
                          res[i],
                          dims,
                          cut_a=cut_a[i], 
                          cut_b=cut_p[i],
                          upper=upper_case[i],
                          w=jw_weight[i],
                          partial=partials[i])
        elseif match_method[i] == "string"
            gammaCKpar!(dfA[!,varnames[i]],
                        dfB[!,varnames[i]],
                        res[i],
                        dims,
                        distmethod=stringdist_method[i],
                        cut_a=cut_a[i], 
                        cut_b=cut_p[i],
                        w=jw_weight[i],
                        partial=partials[i])
        elseif match_method[i] == "exact" || match_method[i] == "bool"
            gammaKpar!(dfA[!,varnames[i]],
                       dfB[!,varnames[i]],
                       res[i],
                       dims)
        elseif match_method == "numeric" || match_method=="float" || match_method == "int"
            gammaNUMCKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           cut_a=cut_a[i],
                           cut_b=cut_p[i],
                           partial=partials[i])
        end
    end
    @info "Getting table counts"
    counts = get_match_patterns(res)

    @info "Running expectation maximization function"
    resultsEM = emlinkMARmov(counts, obs_a,obs_b,
                             varnames,tol=tol_em,
                             prior_lambda=prior_lambda, w_lambda=w_lambda,
                             prior_pi=prior_pi, w_pi=w_pi,
                             address_field=address_field)

    # testing removing uncessessary indices (where no obs exist)
    #remove_no_matched_var_indices(resultsEM)
    # adding uids 
    resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),))
    
    @info "Retrieving matches"
    getMatches(resultsEM,threshold_match=threshold_match)
    
    return (resultsEM) 
end





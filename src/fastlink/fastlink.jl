
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



function create_comparison_function(res,
                                    col::Int,
                                    match_method::String,
                                    stringdist::String,
                                    jw_weight::Float64,
                                    cut_a::T,
                                    cut_p::T,
                                    upper::Bool,
                                    partial::Bool,
                                    fuzzy::Bool,
                                    comparison_level::Int) where T <: Number
    match2 = [true,true]
    match1 = [true,false]
    missingval = [false,true]
    difference_function = -
    if match_method == "string"
        if fuzzy
            match_fun = gammaCKfuzzy!(view(res.result_matrix,:,res.ranges[col]),
                                      res.array_2Dindex,
                                      res.dims,
                                      cut_a=cut_a,
                                      cut_b=cut_p,
                                      upper=upper,
                                      w=jw_weight,
                                      partial=partial,
                                      match2=reshape(match2,(1,2)),
                                      match1=reshape(match1,(1,2)),
                                      missingval=reshape(missingval,(1,2)))
        else
            match_fun = gammaCKpar!(view(res.result_matrix,:,res.ranges[col]),
                                    res.array_2Dindex,
                                    res.dims,
                                    distmethod=stringdist,
                                    cut_a=cut_a,
                                    cut_b=cut_p,
                                    partial=partial,
                                    w=jw_weight,
                                    match2=match2,
                                    match1=match1,
                                    missingval=missingval)
        end       
    elseif match_method == "exact" || match_method == "bool"
        if comparison_level == 1
            match2 = true
        end
        match_fun = gammaKpar!(view(res.result_matrix,:,res.ranges[col]),
                               res.array_2Dindex,
                               res.dims,
                               match2,missingval)
    elseif match_method == "numeric" || match_method=="float" || match_method == "int"
        match_fun = gammaNUMCKpar!(
            view(res.result_matrix,:,res.ranges[col]),
            res.array_2Dindex,
            res.dims,
            cut_a=cut_a,
            cut_b=cut_p,
            partial=partial,
            match2=match2,
            match1=match1,
            missingval=missingval,
            comparison_function=difference_function,
            match_method=match_method)
    end

    return match_fun
end


# Constructor to check types and variable presence in varnames parameter
struct FastLinkVars
    varnames::Vector{String}
    types::Vector{String}
    comparison_funs::Vector{Function}
    function FastLinkVars(varnames::Vector{String},
                          res::ResultMatrix,
                          vartypes::Vector{String},
                          stringdist::Vector{String},
                          jw_weight::Vector{Float64},
                          cut_a::Vector{T},
                          cut_p::Vector{T},
                          upper::Vector{Bool},
                          partial::Vector{Bool},
                          fuzzy::Vector{Bool},
                          comparison_levels::Vector{Int}
                          ) where T <: Number
        comparison_funs=Function[]
        for i in 1:length(comparison_levels)
            push!(comparison_funs,
                  create_comparison_function(res,i,vartypes[i],stringdist[i],jw_weight[i],
                                             cut_a[1],cut_p[i],upper[i],partial[i],fuzzy[i],comparison_levels[i]
                                             ))
        end
        new(varnames, vartypes, comparison_funs)
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
                  varnames::Vector{String};
                  match_method=String[],
                  partials=[true],
                  fuzzy=[false],
                  upper_case=[true],
                  stringdist_method = ["jw"],
                  cut_a = [0.92], cut_p = [0.88],
                  jw_weight = [0.1],
                  tol_em = 1e-05,
                  threshold_match = 0.85,
                  dedupe_matches = true,
                  verbose = false)
    # Allow missing vals in case one has no missing vals

    # dims
    numvars=length(varnames)
    obs_a=nrow(dfA)
    obs_b=nrow(dfB)

    @info "Checking settings for $numvars declared variables."
    partials = check_input_lengths(partials, numvars, "partials")
    upper_case = check_input_lengths(upper_case, numvars, "upper_case")
    fuzzy = check_input_lengths(fuzzy, numvars, "fuzzy")
    jw_weight = check_input_lengths(jw_weight, numvars, "jw_weight")
    cut_a = check_input_lengths(cut_a, numvars, "cut_a")
    cut_p = check_input_lengths(cut_p, numvars, "cut_p")
    stringdist_method = check_input_lengths(stringdist_method, numvars, "stringdist_method")

    
    vartypes, comparison_levels = check_var_types(dfA,dfB,varnames,match_method,partials)
    
    res=ResultMatrix(comparison_levels, (obs_a,obs_b))

    fastlink_settings=FastLinkVars(varnames,res,vartypes,stringdist_method,jw_weight,cut_a,cut_p,upper_case,partials,fuzzy,comparison_levels)

    

    return () -> begin

        # allow missing for comparisons
        allowmissing!(dfA)
        allowmissing!(dfB)


        # iterate through variables and execute function over them
        for i in eachindex(varnames)
            @info "Now matching var $(varnames[i]) using $(match_method[i])"
            fastlink_settings.comparison_funs[i](dfA[!,varnames[i]],dfB[!,varnames[i]])
        end

        @info "Getting table counts"
        counts = tableCounts(view(res.result_matrix,:,:), varnames)

        @info "Running expectation maximization function"
        resultsEM = emlinkMARmov(counts[2], obs_a,obs_b,
                                 varnames,res.ranges,tol=tol_em)

        @info "Retrieving matches"
        matches = getMatches(resultsEM, counts[1], obs_a,threshold_match=threshold_match)
        
        return (resultsEM, matches)
    end
end

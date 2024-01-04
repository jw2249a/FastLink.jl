# Ensure that vartypes match each other. if string then both must be strings
# TODO: implement tighter parameter checks and coercion where needed. 
function check_var_types(x::DataFrame, y::DataFrame, varnames::Vector{String})
    xtypes=eltype.(eachcol(select(x, varnames)))
    ytypes=eltype.(eachcol(select(y, varnames)))
    # checking if precompile will pick up on this
    for (ix,iy,iv) in zip(xtypes,ytypes, varnames)
        if (ix <: Union{Missing,AbstractString}) ⊽ (iy <: Union{Missing,AbstractString})
            throw("*(VAR $iv)*: dfA type $ix does not match dfB type $ix")
        end
    end
    return xtypes
end

# Constructor to check types and variable presence in varnames parameter
struct fastLinkVars
    varnames::Vector{String}
    types::Vector{Type}    
    function fastLinkVars(dfA::DataFrame, dfB::DataFrame,varnames::Vector{String})
        varsA=names(dfA)
        varsB=names(dfB)
        types = []
        ## checking that each var is validly defined
        for i in varnames
            if all(i .!== varsA)
                throw("Missing var $i from dfA")
            end
            if all(i .!== varsB) 
                throw("Missing var $i from dfB")
            end
        end
        vartypes=check_var_types(dfA,dfB,varnames)
        new(varnames, vartypes)
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
                  fuzzy::Bool=false,
                  string_case="upper",
                  stringdist_method = "jw",
                  cut_a = 0.92, cut_p = 0.88,
                  jw_weight = 0.1,
                  tol_em = 1e-05,
                  threshold_match = 0.85,
                  dedupe_matches = true,
                  verbose = false)
    # Allow missing vals in case one has no missing vals
    allowmissing!(dfA)
    allowmissing!(dfB)

    obs_a=nrow(dfA)
    obs_b=nrow(dfB)
    
    vars=fastLinkVars(dfA, dfB, varnames)
    comparison_levels=[2 for i in vars.varnames]
    res=ResultMatrix(comparison_levels, (obs_a,obs_b))
    
    for col in 1:length(vars.varnames)
        if fuzzy
            gammaCKfuzzy!(dfA[!,vars.varnames[col]],
                          dfB[!,vars.varnames[col]],
                          view(res.result_matrix,:,res.ranges[col]),
                          res.array_2Dindex,
                          res.dims,cut_a=cut_a,cut_b=cut_p,
                          upper=string_case == "upper")
        else
            gammaCKpar!(dfA[!,vars.varnames[col]],
                        dfB[!,vars.varnames[col]],
                        view(res.result_matrix,:,res.ranges[col]),
                        res.array_2Dindex,
                        res.dims,cut_a=cut_a,cut_b=cut_p,
                        distmethod=stringdist_method)
        end
    end
    
    counts = tableCounts(view(res.result_matrix,:,:), varnames)
    
    resultsEM = emlinkMARmov(counts[2], obs_a,obs_b,
                             varnames,res.ranges,tol=tol_em)
    
    matches = getMatches(resultsEM, counts[1], obs_a,threshold_match=threshold_match)
    
    return (resultsEM, matches)
end

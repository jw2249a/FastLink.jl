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
- `dfA::DataFrame`: Table of records to be matched.
- `dfB::DataFrame`: Table of records to be matched.
- `varnames::Vector{String}`: A vector that contains the variable names present in both tables.
- `idvar::Tuple{String,String}`: Tuple of unique ids for records in both tables (should differ from each other). Should be formatted as so that the arguments are ordered in ("dfa","dfb") order.
- `term_freq_adjustment::Vector{Bool}`: Determines whether you want the term frequencies for each comparision for a given variable. Note: does not adjust match weight.
- `match_method::Vector{String}`: Specifies the matching method you would like to do. Match methods supported are "string","exact","fuzzy" (jaro-winkler strings only),"numeric","float",and "int")
- `partials::Vector{Bool}` Specifies whether you want to do 2 (true) or 1 (false) comparison levels for a given variable.
- `upper_case::Vector{Bool}` that specifies whether a strings column value is upper or lower (only if `match_method` for column is "fuzzy").
- `stringdist_method::Vector{String}`: String distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming).
- `cut_a::Float`: First lower bound for string distance cutoff.
- `cut_p::Float`: Second lower bound for string distance (if varnames in partial).
- `jw_weight`: Winkler weight for jw string distance.
- `address_field::Vector{Bool}`: Specifies whether a field is an address field.


# Returns
- NamedTuple with these vars
indices        iter_converge  matched_ids    obs_a
obs_b          p_m            p_u            patterns_b
patterns_w     pgamma_jm      pgamma_ju      pgamma_km
pgamma_ku      tf_adj_table   varnames       zeta_j

# Examples
```julia
matched_data = fastLink(dfA, dfB, ["firstname", "lastname", "city"])
"""
function fastLink(dfA::DataFrame, dfB::DataFrame,
                  varnames::Vector{String},
                  idvar::Tuple{String,String};
                  term_freq_adjustment=[false],
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
                  dedupe_matches = true)
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
    term_freq_adjustment = check_input_lengths(term_freq_adjustment, numvars, "term_freq_adjustment")
    stringdist_method = check_input_lengths(stringdist_method, numvars, "stringdist_method")
    
    vartypes, comparison_levels = check_var_types(dfA,dfB,varnames,match_method,partials)

    # results table
    res = [DiBitMatrix(obs_a,obs_b) for _ in varnames]

    # term frequency tables
    tf_table_x = [ones(Float16,dims[1]) for _ in varnames]
    tf_table_y = [ones(Float16,dims[2]) for _ in varnames]
    
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)

    
    # iterate through variables and execute function over them
    for i in eachindex(varnames)
        @info "Now matching var $(varnames[i]) using $(match_method[i])"
        if match_method[i] == "fuzzy"
            if term_freq_adjustment[i]
                gammaCKfuzzy!(dfA[!,varnames[i]],
                              dfB[!,varnames[i]],
                              res[i],
                              dims,
                              view(tf_table_x[i],:),
                              view(tf_table_y[i],:),
                              cut_a=cut_a[i], 
                              cut_b=cut_p[i],
                              upper=upper_case[i],
                              w=jw_weight[i],
                              partial=partials[i])
            else
                gammaCKfuzzy!(dfA[!,varnames[i]],
                              dfB[!,varnames[i]],
                              res[i],
                              dims,
                              cut_a=cut_a[i], 
                              cut_b=cut_p[i],
                              upper=upper_case[i],
                              w=jw_weight[i],
                              partial=partials[i])
            end
        elseif match_method[i] == "string"
            if term_freq_adjustment[i]
                gammaCKpar!(dfA[!,varnames[i]],
                            dfB[!,varnames[i]],
                            res[i],
                            dims,
                            view(tf_table_x[i],:),
                            view(tf_table_y[i],:),
                            distmethod=stringdist_method[i],
                            cut_a=cut_a[i], 
                            cut_b=cut_p[i],
                            w=jw_weight[i],
                            partial=partials[i])
            else
                gammaCKpar!(dfA[!,varnames[i]],
                            dfB[!,varnames[i]],
                            res[i],
                            dims,
                            distmethod=stringdist_method[i],
                            cut_a=cut_a[i], 
                            cut_b=cut_p[i],
                            w=jw_weight[i],
                            partial=partials[i])
            end
        elseif match_method[i] == "exact" || match_method[i] == "bool"
            if term_freq_adjustment[i]
                gammaKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           dims,
                           view(tf_table_x[i],:),
                           view(tf_table_y[i],:))
            else
                gammaKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           dims)
            end
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

    if any(term_freq_adjustment)
        resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),
                                      tf_adj_table = tf_adj_table(resultsEM,varnames,tf_table_x,tf_table_y)))        
        
    else
        resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),))        
    end
    

    
    @info "Retrieving matches"
    getMatches!(resultsEM,threshold_match=threshold_match)
    return (resultsEM) 
end

# version of fast link that i can pass to via named tuples
function fastLink(dfA::DataFrame, dfB::DataFrame;
                  varnames=String[],
                  match_method=String[],
                  idvar=String[],
                  partials=[true],
                  term_freq_adjustment=[false],
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
                  dedupe_matches = true)

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
    term_freq_adjustment = check_input_lengths(term_freq_adjustment, numvars, "term_freq_adjustment")
    stringdist_method = check_input_lengths(stringdist_method, numvars, "stringdist_method")
    
    vartypes, comparison_levels = check_var_types(dfA,dfB,varnames,match_method,partials)

    # results table
    res = [DiBitMatrix(obs_a,obs_b) for _ in varnames]

    # term frequency tables
    tf_table_x = [ones(Float16,dims[1]) for _ in varnames]
    tf_table_y = [ones(Float16,dims[2]) for _ in varnames]
    
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)

    
    for i in eachindex(varnames)
        @info "Now matching var $(varnames[i]) using $(match_method[i])"
        if match_method[i] == "fuzzy"
            if term_freq_adjustment[i]
                gammaCKfuzzy!(dfA[!,varnames[i]],
                              dfB[!,varnames[i]],
                              res[i],
                              dims,
                              view(tf_table_x[i],:),
                              view(tf_table_y[i],:),
                              cut_a=cut_a[i], 
                              cut_b=cut_p[i],
                              upper=upper_case[i],
                              w=jw_weight[i],
                              partial=partials[i])
            else
                gammaCKfuzzy!(dfA[!,varnames[i]],
                              dfB[!,varnames[i]],
                              res[i],
                              dims,
                              cut_a=cut_a[i], 
                              cut_b=cut_p[i],
                              upper=upper_case[i],
                              w=jw_weight[i],
                              partial=partials[i])
            end
        elseif match_method[i] == "string"
            if term_freq_adjustment[i]
                gammaCKpar!(dfA[!,varnames[i]],
                            dfB[!,varnames[i]],
                            res[i],
                            dims,
                            view(tf_table_x[i],:),
                            view(tf_table_y[i],:),
                            distmethod=stringdist_method[i],
                            cut_a=cut_a[i], 
                            cut_b=cut_p[i],
                            w=jw_weight[i],
                            partial=partials[i])
            else
                gammaCKpar!(dfA[!,varnames[i]],
                            dfB[!,varnames[i]],
                            res[i],
                            dims,
                            distmethod=stringdist_method[i],
                            cut_a=cut_a[i], 
                            cut_b=cut_p[i],
                            w=jw_weight[i],
                            partial=partials[i])
            end
        elseif match_method[i] == "exact" || match_method[i] == "bool"
            if term_freq_adjustment[i]
                gammaKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           dims,
                           view(tf_table_x[i],:),
                           view(tf_table_y[i],:))
            else
                gammaKpar!(dfA[!,varnames[i]],
                           dfB[!,varnames[i]],
                           res[i],
                           dims)
            end
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

    if any(term_freq_adjustment)
        resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),
                                      tf_adj_table = tf_adj_table(resultsEM,varnames,tf_table_x,tf_table_y)))        
        
    else
        resultsEM = merge(resultsEM, (matched_ids = indices_to_uids(dfA[!, idvar[1]],dfB[!, idvar[2]],resultsEM.indices),))        
    end
    

    
    @info "Retrieving matches"
    getMatches!(resultsEM,threshold_match=threshold_match)
    return (resultsEM) 
end





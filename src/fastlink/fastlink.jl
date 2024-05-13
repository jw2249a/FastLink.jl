


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
function fastLink(dfA::DataFrame, dfB::DataFrame, config::Dict{String,Any})
    # idvar to Tuple
    idvar=Tuple(config["idvar"])
    # dims
    varnames = retrieve(config,"varname")
    numvars=length(varnames)
    _dims = (nrow(dfA),nrow(dfB))

    # results table
    res = Dict(v=>DiBitMatrix(_dims...) for v in varnames)
    
    # fetch parameters for varnames
    parameters=Dict(v=>fetch_parameters(config,v) for v in varnames)
    
    tf_tables = Dict{String, Vector{Vector{Float16}}}(v=>[ones(Float16,_dims[1]),ones(Float16,_dims[2])] for v in varnames if haskey(parameters[v], "tf_adjust") && parameters[v]["tf_adjust"])
    
    # structure of the expectation maximization function in order of the ability to be executed
    emlink_configuration = parse_configuration(config)
    
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)
    
    for v in varnames
        match_method = lowercase(parameters[v]["method"])
        term_freq_adjustment = retrieve(parameters[v],"tf_adjust") |> x -> !isempty(x) && x |> first
        
        @info "Now matching var $(v) using $(match_method) with tf_adjust: $term_freq_adjustment"
        if term_freq_adjustment
            comparisons_args=namedtuple(remove_keys(parameters[v], ["method", "varname", "tf_adjust", "tf_adjustment_weight"]))
            
            if match_method == "fuzzy"
                gammaCKfuzzy!(dfA[!,v],
                              dfB[!,v],
                              res[v],
                              view(tf_tables[v][1],:),
                              view(tf_tables[v][2],:);
                              comparisons_args...)
            elseif match_method == "string"
                gammaCKpar!(dfA[!,v],
                            dfB[!,v],
                            res[v],
                            view(tf_tables[v][1],:),
                            view(tf_tables[v][2],:);
                            comparisons_args...)
            elseif match_method ∈ keys(STRING_DISTANCE_METHODS)
                gammaCKpar!(dfA[!,v],
                            dfB[!,v],
                            res[v],
                            view(tf_tables[v][1],:),
                            view(tf_tables[v][2],:);
                            distmethod=STRING_DISTANCE_METHODS[match_method],
                            comparisons_args...)
            elseif match_method == "exact" || match_method == "bool"
                gammaKpar!(dfA[!,v],
                           dfB[!,v],
                           res[v],
                           view(tf_tables[v][1],:),
                           view(tf_tables[v][2],:);
                           comparisons_args...)
            elseif match_method == "numeric" || match_method=="float" || match_method == "int"
                gammaNUMCKpar!(dfA[!,v],
                               dfB[!,v],
                               res[v];
                               comparisons_args...)
            end
        else
            comparisons_args=namedtuple(remove_keys(parameters[v], ["method", "varname"]))
            if match_method == "fuzzy"
                gammaCKfuzzy!(dfA[!,v],
                              dfB[!,v],
                              res[v];
                              comparisons_args...)
            elseif match_method == "string"
                gammaCKpar!(dfA[!,v],
                            dfB[!,v],
                            res[v];
                            comparisons_args...)
            elseif match_method ∈ keys(STRING_DISTANCE_METHODS)
                gammaCKpar!(dfA[!,v],
                            dfB[!,v],
                            res[v];
                            distmethod=STRING_DISTANCE_METHODS[match_method],
                            comparisons_args...)
            elseif match_method == "exact" || match_method == "bool"
                gammaKpar!(dfA[!,v],
                           dfB[!,v],
                           res[v];
                           comparisons_args...)
            elseif match_method == "numeric" || match_method=="float" || match_method == "int"
                gammaNUMCKpar!(dfA[!,v],
                               dfB[!,v],
                               res[v];
                               comparisons_args...)
                end
        end
    end     

    results = process_comparisons(res, emlink_configuration, _dims, parameters, tf_tables)

    if length(results)  == 3
        return Dict("ids" => indices_to_uids(dfA[!, config["idvar"][1]],dfB[!, config["idvar"][2]],results[1].indices),
                "resultsEM" => results[2],
                "resultsTF" => results[3])
    else
        return Dict("ids" => indices_to_uids(dfA[!, config["idvar"][1]],dfB[!, config["idvar"][2]],results[1].indices),
                "resultsEM" => results[2])
    end
end


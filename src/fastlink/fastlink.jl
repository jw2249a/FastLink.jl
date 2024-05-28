


"""
Probabilistic record matching using FastLink data-matching algorithm.
Algorithm taken from:
- Enamorado, Ted, Benjamin Fifield, and Kosuke Imai. 2017. fastLink: Fast Probabilistic Record Linkage with Missing Data. Version 0.6.

# Arguments
- `dfA::DataFrame`: Table of records to be matched.
- `dfB::DataFrame`: Table of records to be matched.
- `config::Dict`: Configuration for the match.
# Examples
```julia
matched_data = fastLink(dfA, dfB, config)
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


function fastLink(dfA::DataFrame, dfB::DataFrame, config::Dict{String,Any}, benchmark_N::Bool)
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
    btimes = 
    # allow missing for comparisons
    allowmissing!(dfA)
    allowmissing!(dfB)

    benchtimes = []
    for v in varnames
        starttime = time()
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
        push!(benchtimes, time() - starttime)
    end     
    
    results = process_comparisons(res, emlink_configuration, _dims, parameters, tf_tables)

    if length(results)  == 3
        return Dict("ids" => indices_to_uids(dfA[!, config["idvar"][1]],dfB[!, config["idvar"][2]],results[1].indices),
                "resultsEM" => results[2],
                    "resultsTF" => results[3],
                    "benchtimes" => benchtimes)
    else
        return Dict("ids" => indices_to_uids(dfA[!, config["idvar"][1]],dfB[!, config["idvar"][2]],results[1].indices),
                    "resultsEM" => results[2],
                     "benchtimes" => benchtimes)
    end
end


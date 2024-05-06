module settings

export namedtuple, fetch_parameters, retrieve, parse_configuration, remove_keys

dictkeys(d::Dict) = (collect(Symbol.(keys(d)))...,)
dictvalues(d::Dict) = (collect(values(d))...,)

namedtuple(d::Dict{String,T}) where {T} =
    NamedTuple{dictkeys(d)}(dictvalues(d))

function fetch_parameters(vec::Vector, parameter::String)
    for value in vec
        retvalue = fetch_parameters(value, parameter)
        if !isnothing(retvalue)
            return retvalue
        end
    end
end

function fetch_parameters(dict::Dict, parameter::String)
    for (key, value) in dict
        if key == "varname"
            if value == parameter
                return deepcopy(dict)
            end
        elseif key == "variables" || key == "comparisons"
            return fetch_parameters(value, parameter)
        end
    end
end


# function to retrieve nested keys
function retrieve(vec::Vector, key_of_interest::String)
    returnval = []
    for value in vec
        if value isa Dict
            append!(returnval,retrieve(value, key_of_interest))
        end
    end
    return returnval
end

function retrieve(dict::Dict, key_of_interest::String)
    returnval = []
    for (key, value) in dict
        if key == key_of_interest
            append!(returnval,[value])
        end
        if value isa Dict || value isa Vector
            append!(returnval,retrieve(value, key_of_interest))
        end
    end
    return returnval
end

remove_keys(x::Dict,y) = filter(((c,v),) -> c ∈ collect(setdiff(keys(x), y)), x)

get_config_varname(x::Dict) = haskey(x,"comparisons") ? x["comparisons"]["name"] : x["varname"]

get_calculated_vars(x::Vector) = [i["comparisons"]["name"] for i in x if haskey(i, "comparisons")]

get_config_variables(x::Vector) = [get_config_varname(i) for i in x]

function extract_comparison_info(x::Dict, parent::Union{Nothing,String})
    return Dict("parent"=>parent, "name"=>x["name"],
                "variables" => get_config_variables(x["variables"]),
                "calculated_vars" => get_calculated_vars(x["variables"]),
                "parameters" => namedtuple(remove_keys(x, ["name","variables"])))
end



function parse_configuration(config::Vector, results::Vector{Vector{Dict{String,Any}}}, parent::String, depth::Int)
    append!(results[depth], extract_comparison_info.(config, parent))
    add_child=true
    for i in config
        comparisons=map(x->x["comparisons"], filter(c -> haskey(c, "comparisons"), i["variables"]))
        if length(comparisons) > 0
            if add_child
                push!(results,[])
                add_child=false
            end
            results = parse_configuration(comparisons, results, i["name"],depth+1)
        end
    end
    
    return results
end

function parse_configuration(config::Dict)
    results = [[extract_comparison_info(config["comparisons"], nothing)]]
    push!(results,[])
    comparisons=map(x->x["comparisons"], filter(c -> haskey(c, "comparisons"), config["comparisons"]["variables"]))
    if length(comparisons) > 0
        return [i for i in reverse(parse_configuration(comparisons, results, config["comparisons"]["name"], 2)) if i != []]
    else
        return [i for i in results if i != []]
    end
end

end



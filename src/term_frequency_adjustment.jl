module tf
using DataFrames
import ..settings: retrieve
import ..match1, ..match2

export bf2_to_probability, bf_to_probability, find_tf_pattern_vars, generate_tf_adjustment_dict, get_tf_adjustment_prior_weights, calculate_tf_denom_fuzzy, calculate_tf_denom_exact

# function get_kval(x, i)
#     x_l = length(x)
#     if i == 0
#         i = 1
#     elseif i == 2
#         i = x_l
#     elseif i == 1 && x_l == 3
#         i = length(x) - 1
#     end
#     return x[i]
# end


function generate_pattern_level_structure(N::Int64, varnames::Vector{String}, prior::Float64)
    return Dict((v => zeros(Float64, N)
                 for v in varnames)...,
                "final_weight" => fill(prior, N),
                "ismatch" => falses(N), "final_zetaj" => zeros(Float64, N), "tf_adjusted"=> true)
end

function generate_tf_skeleton(EMOutput::Dict{String,Any}, tf_indices::Vector{Vector{Int64}})
    collect(1:EMOutput["number_of_unique_patterns"]) .|>
        x -> let colindices = tf_indices[x];
            patterns_w = EMOutput["patterns_w"][x,:];
            count = patterns_w.counts;
            varnames = EMOutput["varnames"]
            if isempty(colindices)
                Dict("tf_adjustment_weight" => [0.0],
                     "final_weight" => [Float64(patterns_w.weights)],
                     "final_zetaj" => [bf_to_probability(Float64(patterns_w.weights))],
                     "ismatch" => bf_to_probability(Float64(patterns_w.weights)) >= EMOutput["threshold_match"],
                     "tf_adjusted" => false)
            else
                generate_pattern_level_structure(count, varnames[colindices], Float64(patterns_w.weights))
            end
        end
end

function get_tf_adjustment_prior_weights(parameters::Dict{String,Dict{String,Any}}, tf_vars::Vector{String})
    weights=[retrieve(parameters[v], "tf_adjustment_weight") for v in tf_vars]
    [length(w) == 0 ? 1.0 : w[1] for w in weights]
end


function pattern_tf_adjustment!(tf_result::Dict{String, Any}, colnames::Vector{String},
                                tfvals::Vector{Vector{Float16}}, uvals::Vector{Float64},
                                threshold_match::Float64, count::Int64, prior_weight::Vector{Float64}; base="log2")
    if base == "log2"
        compute_tf_weight = tf_to_bf2
        bf_to_prob = bf2_to_probability
    else
        compute_tf_weight = tf_to_bf
        bf_to_prob = bf_to_probability
    end
    Threads.@threads for ri in collect(1:count)
        for (ci, cn) in enumerate(colnames)
            adjustment = compute_tf_weight(tfvals[ci][ri], uvals[ci], prior_weight[ci])
            tf_result[cn][ri] = adjustment
            tf_result["final_weight"][ri] += adjustment
        end
        zetaj = bf_to_prob(tf_result["final_weight"][ri])
        tf_result["final_zetaj"][ri] = zetaj
        tf_result["ismatch"][ri] = zetaj >= threshold_match
    end
    return nothing
end

function generate_tf_adjustment_dict(EMOutput::Dict{String,Any},e::Dict{String, Any}, tf_vars::Vector{String}, tfPatterns::Dict{String,Vector}, tf_prior_weights::Vector{Float64}; base="log2")
    tfResults = generate_tf_skeleton(EMOutput, tfPatterns["relevant_tf_indices"])
    threshold_match = EMOutput["threshold_match"]
    for pattern_id in collect(1:EMOutput["number_of_unique_patterns"])
        colindices = tfPatterns["relevant_tf_indices"][pattern_id]        
        
        count = EMOutput["patterns_w"].counts[pattern_id]
        tf_uvals = get_tf_u_values(EMOutput["patterns_w"], colindices, pattern_id)

        ci = [findfirst(v .== tf_vars) for v in e["variables"][colindices]]
        tf_pw = tf_prior_weights[ci]
        
        if tfResults[pattern_id]["tf_adjusted"]
            
            pattern_tf_adjustment!(tfResults[pattern_id],
                                   EMOutput["varnames"][colindices],
                                   tfPatterns["tf_denom_vals"][pattern_id],
                                   tf_uvals,
                                   threshold_match,
                                   count,
                                   tf_pw;
                                   base)
            
            
        end
    end
    return tfResults
end


function find_tf_pattern_vars(unique_pattern::Vector{UInt8}, tf_indices::Vector{Int64})::Vector{Int64}
    return intersect(tf_indices, findall(unique_pattern .== match1 .|| unique_pattern .== match2))
end

function calculate_tf_denom_fuzzy(x, y)
    return max(x,y)
end

function calculate_tf_denom_exact(x, y)
    return (x + y)/2
end

function tf_to_bf(denom, uval, prior_weight)
    return log(uval/denom) * prior_weight
end

function tf_to_bf2(denom, uval, prior_weight)
    return log2(uval/denom) * prior_weight
end

function bf_to_probability(w::Float64)::Float64
    return exp(w)/(1+exp(w))
end
function bf2_to_probability(w::Float64)::Float64
    return exp2(w)/(1+exp2(w))
end

function get_tf_u_values(patterns_w::DataFrame, tf_indices::Vector{Int64}, row::Int64)
    N_comparisons = sum(patterns_w[:,"counts"])
    tf_uvals=[[sum(patterns_w[ismissing.(patterns_w[:,tfi]) .== false .&& patterns_w[:,tfi] .== value,"counts"])/N_comparisons for value in 1:2] for tfi in tf_indices]
    return [tf_uvals[i][v] for (i,v) in enumerate(patterns_w[row,tf_indices])]
end

# applies term frequency adjustments to table
function tf_adj_table(resultsEM::NamedTuple,varnames::Vector{String},tf_table_x::Vector{Vector{Float16}},tf_table_y::Vector{Vector{Float16}})
    tf_vec = [DataFrame() for _ in eachindex(resultsEM.indices)]
    new_names=vcat("tf_" .* varnames .* "_a", "tf_" .*  varnames .* "_b")
    for i in eachindex(resultsEM.indices)
        result_len=length(resultsEM.indices[i])
        tf_results=DataFrame(ones(Float16,(result_len, 2*length(varnames))),new_names)
        Threads.@threads for ii in 1:result_len
            val=resultsEM.indices[i][ii]
            rowval=vcat([tf_table_x[varid][val.row] for varid in eachindex(varnames)],[tf_table_y[varid][val.col] for varid in eachindex(varnames)])
            tf_results[ii,:] = rowval
        end

        tf_vec[i] = tf_results
    end

    return tf_vec
end

function get_tf_weight(parameters::Dict, varname::String)::Float64
    return in("tf_adjustment_weight", keys(parameters[varname])) ? parameters[varname]["tf_adjustment_weight"] : 1.0
end


function update_tf_zetas_ismatch(final_bf::Vector{Float64},
                                 threshold_match::Float64)
    bf_len = length(final_bf)
    final_zj = zeros(Float64, bf_len)
    ismatch = falses(bf_len)
    Threads.@threads for (i, v) in collect(zip(1:bf_len,final_bf))
        final_zj[i] = bf_to_probability(v)
        ismatch[i] = final_zj[i] > threshold_match
    end
    return (final_zj, ismatch)
end

end

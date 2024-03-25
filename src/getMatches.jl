# Decomposes match vector index back into original indices from initial inputs.
function get_2Dindex(index, nrows)
    zero_based_index = index - 1
    row = Int(mod(zero_based_index, nrows)) + 1
    col = Int(div(zero_based_index, nrows)) + 1
    return (row, col)
end

"""
Converts the matches from the tableCounts function based on the predefined threshold from the fastLink call into original indices from the input datasets. 

# Arguments
- `resultsEM::NamedTuple`: Output of the expectation maximization fuction (eg emlinkMARmov())
- `threshold_match`: Lower bound for the posterior probability that will act as a cutoff for matches.
"""
function getMatches!(resultsEM::NamedTuple;
                    threshold_match=0.85,u_b=1e10)
    resultsEM.patterns_w.ismatch = resultsEM.zeta_j .>= threshold_match .&& resultsEM.patterns_w.weights .<= u_b
    return nothing
end

# applies term frequency adjustments to table
function tf_adj_table(resultsEM::NamedTuple,varnames::Vector{String},tf_table_x::Vector{Vector{Float16}},tf_table_y::Vector{Vector{Float16}})
    tf_vec = [DataFrame() for _ in eachindex(resultsEM.indices)]
    new_names=vcat("tf_" .* varnames .* "_x", "tf_" .*  varnames .* "_y")
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

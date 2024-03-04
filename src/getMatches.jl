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
function getMatches(resultsEM::NamedTuple;
                    threshold_match=0.85,u_b=1e10)
    resultsEM.patterns_w.ismatch = resultsEM.zeta_j .>= threshold_match .&& resultsEM.patterns_w.weights .<= u_b
    return nothing
end



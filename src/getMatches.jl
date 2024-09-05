module getmatches
using DataFrames
"""
Function to extract matches from fastlink output. Specify an threshold for matches or use the default threshold specified in the configuration for the fastlink function. returns a `Vector{DataFrame}` that can be matched to the same indices in the patterns_w object subset by zeta_j >= threshold_match.
"""
function getMatches(FastLinkOutput::Dict{String, Any}, threshold_match::T) where T <: AbstractFloat
    (FastLinkOutput["resultsEM"]["patterns_w"].zeta_j .>= threshold_match) |>
        findall |>
        ids -> FastLinkOutput["ids"][ids] .|>
        x->DataFrame(x,FastLinkOutput["idvar"])
end
function getMatches(FastLinkOutput::Dict{String, Any})
    getMatches(FastLinkOutput,FastLinkOutput["resultsEM"]["threshold_match"])
end

export getMatches
end # module getMatches

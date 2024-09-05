# FastLink.jl
Fast Probabilistic Record Linkage for the Julia Language
## What is FastLink.jl

The purpose of FastLink.jl is to bring a fast record linkage package to the julia language. When attempting to match large datasets using existing libraries in R and Python, I found they can be very slow and succumb to issues with memory pressure. This implementation of the fastlink algorithm is intended to scale effeciently in parallel and be able to easily handle matches between tabular data that span millions of rows. 

[![Run tests](https://github.com/jw2249a/FastLink.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jw2249a/FastLink.jl/actions/workflows/test.yml)

___________________________
### Using the fastLink function

The basic arguments for the `fastLink` function to run are

- `dfA`: A `DataFrame` table of records to be matched.

- `dfB`: A `DataFrame` table of records to be matched.

- `config`: A `Dict{String, Any}` that specifies how the two dataframes ought to be matched. 

### Match Configurations

The match configuration for a FastLink match needs to contain certain in the base dictionary (nested dictionaries will be discussed later).

The Base Dictionary needs to contain:

- `link_type`: Either `link_only`, `dedupe_only`, or `link_and_dedupe`.

- `idvar`: A `Vector{String}` of length 2 that specifies the ids of the two dataframes (in `[dfA, dfB]` order).

- `comparisons`: a `Dict{String, Any}` a that defines the type of matching to be done and the variables that will be matched. 

### Comparisons dictionary
The comparison dictionary defined above can be located in the base Dictionary or can be substituted instead of a `varname` dictionary in the `variables` vector. The effect of nesting the comparisons in the `variables` vector will lead it to be matched first using the fastlink algorithm and then treated as a single variables in the parent `comparisons` dictionary. You can substitute multiple `varnames` for comparisons at the same level of nestedness. 

Each `comparisons` dictionary much have: 
- `name`: should be "total" in the base dictionary and then can be any `name` for the nested dictionaries. 

- `variables`: a `Vector{Dict{String, Any}}` that contains the individual variable dictionaries and/or a `comparisons` dictionaries.

The optional parameters for the `comparisons` dictionary are:
- `w_lambda::Float64`: Default 0.0.

- `prior_lambda::Float64`: Default 0.0.

- `threshold_match`: Lower bound for the posterior probability that will act as a cutoff for matches. Default `[0.85]`.

- `tol_em`: Convergence tolerance for the EM Algorithm. (default `1e-05`)

- `prior_pi::Float64`: Default 0.0.

- `w_pi::Float64`: Default 0.0.

### Variables

Individual variables can be declared in a dictionary and must contain both a `varname` and `method`. 

- `varname`: name of the variable in `dfA` and `dfB` to be compared.

- `method`: the method to match the variable. The current accepted methods are (`exact`, `fuzzy`, `string`, `numeric`, `float`, `int` any of the `distmethod` options).

### Methods
Each `method` has a number of arguments that can be specified for that matching method. 

#### fuzzy

#### string

#### numeric


- `term_freq_adjustment`: A `Bool` that determines whether you want the term frequencies for each comparision for a given variable. Note: does not adjust match weight.

- `tf_adjustment_weight`: how much to weight on the term_freq_adjustment vs the predicted match value.

- `tf_minimum_u_value`: minimum term frequency value to adjust by.

- `partial`: A `Bool` that specifies whether you want to do 2 (true) or 1 (false) comparison levels for a given variable. Default value `true`. 

- `upper_case`: A `Bool` that specifies whether a strings column value is upper or lower (only if `method`=`true`. Default value is `true`.

- `stringdist_method`: A `String` that specifies the desired string distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming). Default `"jw"`.

- `cut_a`:  A `Float` that specifies the first lower bound for string distance cutoff for each comparison. Default `0.92`.
- `cut_b`: A `Float` that specifies the second lower bound for string distance (if varnames in partial) for each comparison. Default `0.88`.

- `w`: A `Float` that specifies the Winkler weight for jw string distance for each comparison. Default `0.1`.

## Example FastLink configuration with 1 embedded comparisons dictionary 

``` json
{
    "link_type": "link_only",
    "idvar": ["id", "id2"],
    "comparisons": {
        "name": "total",
        "prior_lambda": 0.000001,
        "w_lambda": 0.5,
        "threshold_match": 0.88,
        "variables": [
            {"varname": "firstname", "method": "fuzzy", "partial": true, "cut_a": 0.92, "cut_b": 0.88, "upper": true, "tf_adjust": true, "w": 0.1},
            {"varname": "middlename", "method": "exact"},
            {"varname": "lastname", "method": "jarowinkler", "tf_adjust": true},
            {"varname": "birthyear", "method": "exact"},
            {
                "comparisons": {
                    "name": "address",
                    "threshold_match": 0.92,
                    "variables": [
                        {"varname": "housenum", "method": "exact", "tf_adjust": true},
                        {"varname": "streetname", "method": "jarowinkler", "w": 0.1, "tf_adjust": true, "tf_adjustment_weight":0.25, "tf_minimum_u_value": 0.001},
                        {"varname": "city", "method": "jarowinkler", "tf_adjustment_weight":0.15, "tf_adjust": true}
                    ]
                }
            }
        ]
    }
}
```

__________________
### `fastLink`'s output

For ease of extracting matches, the `getMatches` function was added. You can call it on the fastLink output as the single argument `getMatches(FastLinkOutput)` or with a specified threshold `getMatches(FastLinkOutput, threshold_match)`.

The FastLink output is:
A `Dict{String,Any}` with these vars:
- `ids`: A vector of vectors of tuple pairs of ids for each match pattern.
- `idvar`: ID variable from configuration
- `resultsEM`: The results of the Expectation Maximization algorithm

If term frequency is specified then 
- `resultsTF`: term frequencies for each variable with specified term frequency by pattern if relevant for the pattern (if no term frequency is applied then tf_adjusted is false).

If benchmark is specified:
- `benchtimes`: times for each variable to be matched.



Within `resultsEM` in the EM output, there is:
- `iter_converge` - number of iterations for expectation maximization algorithm to converge. 

- `obs_a` - observations in `dfA`

- `obs_b` - observations in `dfB`

- `p_m` - posterior match probability

- `p_u` - posterior **not** match probability

- `number_of_unique_patterns` - equivalent to number of rows in `patterns_w`

- `number_of_comparisons` - For convenience `nrow(dfA) * nrow(dfB)`

- `patterns_w` - a `DataFrame` of:
  - `gamma_*` - An `Int64` with the gamma values for each variable (similar to `patterns_b`)
  - `counts` - An `Int64` with counts for each agreement pattern
  - `weights` - An `Int64` with partial match weights for each agreement pattern
  - `p_gamma_jm` - A `Float64` that has the posterior probability that a pair matches for each agreement pattern
  - `p_gamma_ju` - A `Float64` that has the posterior probability that a pair **does not** match for each agreement pattern
  - `is_match` - A `Bool` that specifies whether the given pattern is above the input parameter `threshold_match`

- `patterns_b` - vector of all patterns observed. each pattern as a scored number for each variable (0 nonmatch, 1 partial, 2 exact, 3 missing)

- `pgamma_km` - A `Vector{Vector{Float64}}` with posterior probababilities for each variable in the EM algorithm. Ordered (0,1,2).

- `pgamma_ku` - A `Vector{Vector{Float64}}` with posterior probababilities for each variable in the EM algorithm. Ordered (2,1,0).

- `p_gamma_jm` - A `Float64` that has the posterior probability that a pair matches for each agreement pattern  (see `patterns_w`).
- `p_gamma_ju` - A `Float64` that has the posterior probability that a pair **does not** match for each agreement pattern (see `patterns_w`).

- `varnames` - A `Vector{String}` of the input variable names


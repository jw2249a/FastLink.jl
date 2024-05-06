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

- `varnames` A `Vector{String}` that contains the variable names present in both tables.

- `idvar`: A `Tuple{String,String}` that has the unique ids for records in both tables (should differ from each other). Should be formatted as so that the arguments are ordered in ("dfa","dfb") order.

In addition, there are a number of optional parameters to assist with record linkage. For each of the parameters, you can either specify each variables value or specify a vector with 1 value to apply it to all relevant variables. If irrelevant, like `stringdist_method` for a numeric match method, it will be ignored. 

The optional parameters are

- `term_freq_adjustment`: A `Vector{Bool}` that determines whether you want the term frequencies for each comparision for a given variable. Note: does not adjust match weight. Default value `[false]`.

- `match_method`: A `Vector{String}` that specifies the matching method you would like to do. Match methods supported are "string","exact","fuzzy" (jaro-winkler strings only),"numeric","float",and "int"). Default value inferred from column type.

- `partials`: A `Vector{Bool}` that specifies whether you want to do 2 (true) or 1 (false) comparison levels for a given variable. Default value `[true]`. 

- `upper_case`: A `Vector{Bool}` that specifies whether a strings column value is upper or lower (only if `match_method` for column is "fuzzy"). Default value is `[true]`.

- `stringdist_method`: A `Vector{String}` that specifies the desired string distance method ("jw" Jaro-Winkler (Default), "dl" Damerau-Levenshtein, "jaro" Jaro, "lv" Levenshtein, and "ham" Hamming). Default `["jw"]`.

- `cut_a`  A `Vector{Float}` that specifies the first lower bound for string distance cutoff for each comparison. Default `[0.92]`.

- `cut_p`: A `Vector{Float}` that specifies the second lower bound for string distance (if varnames in partial) for each comparison. Default `[0.88]`.

- `jw_weight`: A `Vector{Float}` that specifies the Winkler weight for jw string distance for each comparison. Default `[0.1]`.

- `address_field`: A `Vector{Bool}` that specifies whether a comparison contains an address field. Default `[false]`.

- `tol_em`: Convergence tolerance for the EM Algorithm. (default `1e-05`)

- `threshold_match`: Lower bound for the posterior probability that will act as a cutoff for matches. Default `[0.85]`.

- `prior_lambda::Float64`: Default 0.0.

- `w_lambda::Float64`: Default 0.0.

- `prior_pi::Float64`: Default 0.0.

- `w_pi::Float64`: Default 0.0.

- `dedupe_matches`: Whether to dedupe the matches within the dataset.


__________________
### `fastLink`'s output

A `NamedTuple` with these vars:

- `indices` - a vector with indices in `dfA` and `dfB` that are in each pattern group (see `patterns_w` or `patterns_b`)

- `matched_ids` - same as `indices` but using `idvars` from input parameters

- `iter_converge` - number of iterations for expectation maximization algorithm to converge. 

- `obs_a` - observations in `dfA`

- `obs_b` - observations in `dfB`

- `p_m` - posterior match probability

- `p_u` - posterior **not** match probability

- `patterns_w` - a `DataFrame` of:
  - `gamma_` - An `Int64` with the gamma values for each variable (similar to `patterns_b`)
  - `counts` - An `Int64` with counts for each agreement pattern
  - `weights` - An `Int64` with partial match weights for each agreement pattern
  - `p_gamma_jm` - A `Float64` that has the posterior probability that a pair matches for each agreement pattern
  - `p_gamma_ju` - A `Float64` that has the posterior probability that a pair **does not** match for each agreement pattern
  - `is_match` - A `Bool` that specifies whether the given pattern is above the input parameter `threshold_match`

- `patterns_b` - vector of all patterns observed. each pattern as a scored number for each variable (0 nonmatch, 1 partial, 2 exact, 3 missing)

- `pgamma_km` - A `Vector{Vector{Float64}}` with posterior probababilities for each variable in the EM algorithm. Ordered (0,1,2).

- `pgamma_ku` - A `Vector{Vector{Float64}}` with posterior probababilities for each variable in the EM algorithm. Ordered (2,1,0).

- `tf_adj_table` - A `Vector{DataFrame}` that has a DataFrame for each match pattern and a row in each DataFrame for each comparison appended with the letter of their corresponding dataset.

- `varnames` - A `Vector{String}` of the input variable names
 
- `zeta_j` - A `Vector{Float64}` with the posterior match probabilities for each agreement pattern. 

# Examples
```julia
matched_data = fastLink(dfA, dfB, ["firstname", "lastname", "city"], ("id","id2"))
``` 

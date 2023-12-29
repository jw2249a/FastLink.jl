function fastLink(dfA::DataFrame, dfB::DataFrame,
                  varnames::Vector{String};
                  stringdist_method = "jw",
                  cut_a = 0.94, cut_p = 0.88,
                  jw_weight = 0.1,
                  cut_a_num = 1, cut_p_num = 2.5,
                  tol_em = 1e-04,
                  threshold_match = 0.85,
                  return_all = false,
                  cond_indep = true,
                  estimate_only = false,
                  dedupe_matches = true,
                  return_df = false,
                  linprog_dedupe = false,
                  verbose = false)


    function check_for_variables(dfA, dfB, varnames)

    end
    
    if !all(map(x-> any(x .== names(dfA)), varnames))
        throw("")
    end
    
    result_matrix
    for var in varnames

    end


end

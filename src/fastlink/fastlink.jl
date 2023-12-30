function check_var_types(x::DataFrame, y::DataFrame, varnames::Vector{String})
    xtypes=eltype.(eachcol(select(x, varnames)))
    ytypes=eltype.(eachcol(select(y, varnames)))

    for (ix,iy,iv) in zip(xtypes,ytypes, varnames)
        if ix !== iy
            throw("*(VAR $iv)*: dfA type $ix does not match dfB type $ix")
        end
    end
    return xtypes
end


struct fastLinkVars
    varsnames::Vector{String}
    types::Vector{Type}
    
    function fastLinkVars(dfA::DataFrame, dfB::DataFrame,varnames::Vector{String})

        varsA=names(dfA)
        varsB=names(dfB)
        types = []
        ## checking that each var is validly defined
        for i in varnames
            if all(i .!== varsA)
                throw("Missing var $i from dfA")
            end
            if all(i .!== varsB) 
                throw("Missing var $i from dfB")

            end
        end
        vartypes=check_var_types(dfA,dfB,varnames)
        new(varnames, vartypes)
    end  
end

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
    # Allow missing vals in case one has no missing vals
    allowmissing!(dfA)
    allowmissing!(dfB)
    
    vars=fastLinkVars(dfA, dfB, varnames)
    comparison_levels=[2 for i in varnames]
    res=ResultMatrix(comparison_levels, (nrow(dfA), nrow(dfB)))
    for col in 1:length(varnames)
        gammaCKpar!(dfA[!,varnames[col]],
                    dfB[!,varnames[col]],
                    view(res.result_matrix,:,res.ranges[col]),
                    res.array_2Dindex,
                    res.dims)
    end

    counts = tableCounts(res.result_matrix)

    return counts

end

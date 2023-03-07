cormat_nf = function(df){
    ## Define correlation matrix and pvalue matrix
    rhomat = pvalmat = matrix(NA, ncol(df), ncol(df))

    ## Make the diagonals = 1
    diag(rhomat) = diag(pvalmat) = 1
    
    ## For Loop
    for (i in 1:ncol(df)) {
        for (j in 1:ncol(df)) {
            x = unlist(df[, i])
            y = unlist(df[, j])
            ## define a temporary data frame
            temp_df = data.frame(x = x, y = y)
            ## remove missingness
            temp_df = na.omit(temp_df)
            ## test for data size
            n = nrow(temp_df)
            ## test for variation
            v = apply(temp_df,2,function(x){ length(unique(x)) } ) 
            v = sum( v >=2 )
            ## IF the data set has at least 20 observations
            ##    AND there is variation in both columns
            ##    continue
            if(n >= 20 & v == 2){
              if (class(x) == "factor" & class(y) == "factor") {
                test = factor_on_factor(x, y)
                rhomat[i, j] = test[1]
                pvalmat[i, j] = test[3]
                }
                else {
                    if (class(x) == "numeric" & class(y) == "numeric") {
                      test = numeric_on_numeric(x, y)
                      rhomat[i, j] = test[1]
                      pvalmat[i, j] = test[2]
                    }
                    else {
                      if (class(x) == "factor") {
                        test = factor_on_numeric(cat_values = x, 
                          num_values = y)
                        rhomat[i, j] = test[2]
                        pvalmat[i, j] = test[4]
                      }
                      else {
                        test = factor_on_numeric(cat_values = y, 
                          num_values = x)
                        rhomat[i, j] = test[2]
                        pvalmat[i, j] = test[4]
                      }
                    }
                  }
                } else {
                  rhomat[i, j] = 0
                  pvalmat[i, j] = 1
                }
              }
            }

            
    out = list(RhoMat = rhomat, PvalMat = pvalmat)
    return(out)
}
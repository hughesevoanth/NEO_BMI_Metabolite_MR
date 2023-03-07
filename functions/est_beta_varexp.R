est_beta_varexp = function( wdata,
                  response_variable = NA,
                  covariates = NA,
                  instrument = NA,
                  weights = NA){

  ######################
  ### I. Define Univariate 
  ###      Linear model
  ######################
  univariate_form = formula(paste0(response_variable, " ~ ", instrument ))
  
  ######################
  ### II. Define Multivariate 
  ###      Linear model
  ######################
  if( length(na.omit(covariates)) > 0 ){
    multivariate_form = formula(paste0(response_variable, " ~ ", paste0( covariates, collapse = " + ") , " + ", instrument ))
  } else {
    multivariate_form = NA
  }
    
  #########################
  ## III. RUN Generalized
  ##      Univariate LINEAR MODEL
  #########################
  if( is.na(weights) ){
    uni_mod = glm(univariate_form, data = wdata, family = gaussian() )
  } else {
    ## re-name weights column
    w = grep( weights, colnames(wdata) ); colnames(wdata)[w] = "study_weights"
    uni_mod = glm(univariate_form, weights = study_weights, data = wdata,  family = gaussian() )  
  }
  
  #########################
  ## IV. RUN Generalized
  ##      Multivariate LINEAR MODEL
  #########################
  if( length(multivariate_form)>1 ){
    if( is.na(weights) ){
      multi_mod = glm(multivariate_form, data = wdata, family = gaussian() )
    } else {
      ## weights re-named above in step IV
      multi_mod = glm(multivariate_form, weights = study_weights, data = wdata, family = gaussian() )  
    }
  } else {
    multi_mod = NULL
  }
  
  #########################
  ## V. Beta & Deviance Explained
  #########################
  uni_res = residuals(uni_mod)
  uni_W = normW( uni_res ); names(uni_W) = "uni_W"
  ## Univariate Estimate
  uni_coefs = summary(uni_mod)$coef[instrument,c(1,2,4)]
  names(uni_coefs) = c("uni_beta","uni_se","uni_P")
  ## Dev explained
  a = anova(uni_mod)
  dev = c(a[2,2], a[2,4])
  names(dev) = c(rownames(a)[-1], "residual")
  eta_sq = dev / sum(dev)
  univ_devexp = eta_sq[instrument]
  names(univ_devexp) = "univ_dev_exp"
  
  ## Multivariate Estimate
  if( length(multi_mod) > 0 ){
    ## normality of residuals
    multi_res = residuals(multi_mod)
    multi_W = normW( multi_res ); names(multi_W) = "multi_W"
    ## Univariate Estimate
    multi_coefs = summary(multi_mod)$coef[instrument,c(1,2,4)]
    ##
    a = anova(multi_mod)
    dev = c(a[-1,2], a[nrow(a),4] )
    names(dev) = c(rownames(a)[-1], "residual")
    eta_sq = dev / sum(dev)
    multi_devexp = eta_sq[instrument]
  } else {
    multi_W = NA
    multi_coefs = c(NA, NA, NA)
    multi_devexp = NA
  }
  names(multi_W) = "multi_W"
  names(multi_devexp) = "multi_dev_exp"
  names(multi_coefs) = c("multi_beta","multi_se","multi_P")
  
  
  #########################
  ## VI. Return Data
  #########################
  names(response_variable) = "response_variable"
  names(instrument) = "instrument"
  
  out = c(response_variable, instrument,  
          uni_W, uni_coefs, univ_devexp, 
          multi_W, multi_coefs, multi_devexp
          ) 

  return(out)

}

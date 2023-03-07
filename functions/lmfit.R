lmfit = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  covariates = NA,
                  weights = NA,
                  rnt_outcome = FALSE){

  ############################################
  ## I. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    wdata[, outcome] = rntransform( wdata[, outcome] )
  }

  ######################
  ### II. Define Linear model
  ######################
  if( length(na.omit(covariates)) > 0 ){
    form = formula(paste0(outcome, " ~ ", paste0( covariates, collapse = " + ") , " + ", exposure ))
  } else {
    form = formula(paste0(outcome, " ~ ", exposure ))
  }

  #########################
  ## III. RUN Generalized LINEAR MODEL
  #########################
  if( is.na(weights) ){
    glm_mod = glm(form, data = wdata, family = gaussian() )
  } else {
    ## re-name weights column
    w = grep( weights, colnames(wdata) ); colnames(wdata)[w] = "study_weights"
    glm_mod = glm(form, weights = study_weights, data = wdata,  family = gaussian() )  
  }
  
  #########################
  ## IV. Summary Stats
  #########################
  ## sample size in model
  res = residuals(glm_mod)
  n = length(res)
  names(n) = "n"

  ## normality of residuals
  Wstat = normW(res)
  names(Wstat) = "W"
  
  ## Breusch-Pagan test of homoskedasticity
  Breusch_Pagan_P = lmtest::bptest(glm_mod)$p.value
  names(Breusch_Pagan_P) = "Breusch_Pagan_P"
  
  ## model summary
  s = summary(glm_mod)
  
  ## model coefficients
  glm_coef = s$coefficients
  
  ## Report beta, se, t-value, and P-value
  glm_estimates = glm_coef[exposure, ]; names(glm_estimates) = c("beta","se","tval","P")
  
  ## extract deviances
  a = anova(glm_mod)
  deviance = c( na.omit(a[,2]), a[nrow(a),4] )
  names(deviance) = c(rownames(a)[-1], "residual")
  
  ## estimate eta-sq or deviance explaine
  eta_sq = deviance/sum(deviance)
  
  devexp_by_model = 1 - eta_sq[length(eta_sq)]
  names(devexp_by_model) = "devexp_by_model"
  
  devexp_by_exposure = eta_sq[exposure]
  names(devexp_by_exposure) = "devexp_by_exposure"
  
  ## Sandwich co-variance matrix
  ## ** removes the assumption that residual errors have constant variance **
  vcov_mat = sandwich::vcovHC(glm_mod, type = "HC") ## HC = White’s estimator
  ## Sandwich SE's
  sandwich_se <- diag(vcov_mat)^0.5
  ## re-estimate p_values and make Coefficient table
  coef_table = lmtest::coeftest(glm_mod, vcov = vcov_mat)
  sandwich_estimates = coef_table[exposure, ]
  names(sandwich_estimates) = paste0("sw_", c("beta","se","zval","P"))
  
  ######################
  ## V. Linear model data out
  ######################
  names(exposure) = "exposure"
  names(outcome) = "outcome"
  
  glm_out = c(exposure, outcome, n, Wstat, Breusch_Pagan_P, devexp_by_model, devexp_by_exposure, glm_estimates, sandwich_estimates ) 

  return(glm_out)

}

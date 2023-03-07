ivregfit = function( wdata,
                     outcome,
                     exposure,
                     instrument,
                     covariates = NA,
                     weights = NA,
                     rnt_outcome = FALSE){

  ############################################
  ## I. rank normalize the dependent|outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    wdata[, outcome] = rntransform( wdata[, outcome] )
  }

  ######################
  ### II. Build model
  ######################
  if( length(na.omit(covariates))>0 ){
    form = formula(paste0( outcome, " ~ ", paste0( covariates, collapse = " + ") ," | ", exposure, " | " , instrument ))
  } else {
    form = formula( paste0( outcome, " ~ ",  exposure, " | " , instrument ) )
  }

  #########################
  ## III. RUN MR model
  #########################
  if( is.na( weights) ){
    iv_mod = ivreg(form, data = wdata )
  } else {
    # iv_mod = ivreg( form, data = wdata, weights = wdata[ , weights] )
    w = grep( weights, colnames(wdata) ); colnames(wdata)[w] = "study_weights"
    iv_mod = ivreg( form, weights = study_weights, data = wdata )
  }
  
  ######################
  ### IV. Summary Stats
  ######################
  ## sample size in model
  res = residuals(iv_mod)
  n = length(res); names(n) = "n"
  
  ## normality of residuals
  Wstat = normW(res)
  names(Wstat) = "Wstat"
  
  ## Breusch-Pagan test of homoskedasticity
  Breusch_Pagan_P = lmtest::bptest(iv_mod)$p.value
  names(Breusch_Pagan_P) = "Breusch_Pagan_P"
  
  ## IV summary
  s = summary(iv_mod)
  ## model coefficients
  coef = s$coefficients
  ## Report beta, se, t-value, and P-value
  estimates = coef[exposure,]; names(estimates) = c("beta","se","tval","P")
  ## R-squared
  rsq = s$r.squared
  names(rsq) = "Rsq"
  ## Wald Test for
  wald = s$waldtest[c(1,3,4,2)]
  names(wald) = c("Wald_stat","Wald_df1","Wald_df2","Wald_P")
  ## Diagnostic Test
  diag_test = s$diagnostics
  ## Weak Instrument (F) Test
  WeakInst = diag_test[1, c(3,1,2,4)]
  names(WeakInst) = c("WeakInst_Fstat","WeakInst_df1","WeakInst_df2","WeakInst_P")
  ## Wu Hausman Endogeneity Test
  Wu_Hausman = diag_test[2, c(3,1,2,4)]
  names(Wu_Hausman) = c("Wu_Hausman_stat","Wu_Hausman_df1","Wu_Hausman_df2","Wu_Hausman_P")

  ## Sandwich co-variance matrix
  ## ** removes the assumption that residual errors have constant variance **
  vcov_mat = sandwich::vcovHC(iv_mod, type = "HC") ## HC = White’s estimator
  ## Sandwich SE's
  sandwich_se <- diag(vcov_mat)^0.5
  ## re-estimate p_values and make Coefficient table
  coef_table = lmtest::coeftest(iv_mod, vcov = vcov_mat)
  sandwich_estimates = coef_table[exposure, ]
  names(sandwich_estimates) = paste0("sw_", c("beta","se","zval","P"))
  
  
  ######################
  ## V. MR sum stats out
  ######################
  names(exposure) = "exposure"
  names(outcome) = "outcome"
  
  out = c(exposure, outcome, n, Wstat, Breusch_Pagan_P, rsq, wald, WeakInst, Wu_Hausman, estimates, sandwich_estimates)
  
  return(out)

}

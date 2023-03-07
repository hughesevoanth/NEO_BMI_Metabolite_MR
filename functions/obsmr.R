obsmr = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  instrument = NA,
                  covariates = NA,
                  weights = NA,
                  rnt_outcome = FALSE,
                  outlier_method = "iqr",
                  outlier_cutoff = 5,
                  messages = FALSE){
  
  ############################################
  ### 0. look for any errors in parameters
  ############################################
  # if (!strata[1] %in% c("quartiles","deciles") & length(strata) < 3 & class(strata) != "numeric" ){
  #   stop("strata parameter can either be defined as (1) quartiles, (2) deciles, or (3) a numeric vector of at least length 3 defining strata boundries")
  #  }
  
  ############################################
  ### I. Define Model Data
  ############################################
  # names(outcome) = "outcome"
  if(messages == TRUE){ message("1. defining model data frame.") }
  
  model_variables = na.omit( c(outcome, exposure, instrument, covariates, weights) )
  mod_data = wdata[, c(model_variables)]
  
  ############################################
  ## II. Identify outcome outliers
  ############################################
  if(messages == TRUE){ message("2. identifying outliers") }
  outliers = id_outliers( y = mod_data[, outcome], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)
  
  # How many outcome outliers
  number_of_outcome_outliers = length(outliers); names(number_of_outcome_outliers) = "number_of_outcome_outliers"
  
  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, outcome] = NA
  }
  
  ############################################
  ## III. Identify exposure outliers
  ############################################
  outliers = id_outliers( y = mod_data[, exposure], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)
  
  # How many exposure outliers
  number_of_exposure_outliers = length(outliers); names(number_of_exposure_outliers) = "number_of_exposure_outliers"
  
  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, exposure] = NA
  }
  
  ############################################
  ## IV. Remove NA rows from data set
  ############################################
  mod_data = na.omit(mod_data)
  sample_size = nrow(mod_data)
  names(sample_size) = "sample_size"
  
  ############################################
  ## IV. normality of outcome
  ############################################
  if(messages == TRUE){ message("3. estimating Shapiro-Wilk normality W-stat") }
  W_outcome = normW(mod_data[, outcome]); names(W_outcome) = "W_outcome"
  
  ############################################
  ## IV. normality of exposure
  ############################################
  W_exposure = normW(mod_data[, exposure]); names(W_exposure) = "W_exposure"
  
  ############################################
  ## V. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    if(messages == TRUE){ message("4. Performing rank normal transformation of outcome") }
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  } else {
    if(messages == TRUE){ message("4. No rank normal transformation of outcome performed") }
  }
  
  ######################
  ### VI. Linear model
  ######################
  if(messages == TRUE){ message("5. running observational generalized linear model") }
  obs_stats = lmfit( wdata = mod_data,
                  outcome = outcome,
                  exposure = exposure,
                  covariates = covariates,
                  weights = weights,
                  rnt_outcome = FALSE)
  ## add some name specificity to these stats
  names(obs_stats) = paste0("Obs_", names(obs_stats))
  
  ######################
  ### VII. IVreg model
  ######################
  if(messages == TRUE){ message("6. running TSLS with ivreg") }
  mr_stats = ivregfit( wdata = mod_data,
                      outcome = outcome,
                      exposure = exposure,
                      instrument = instrument,
                      covariates = covariates,
                      weights = weights,
                      rnt_outcome = FALSE)
  ## add some name specificity to these stats
  names(mr_stats) = paste0("MR_", names(mr_stats))
  
  ######################
  ### VIII. IV on exposure 
  ###      beta & variance explained
  ######################
  if(messages == TRUE){ message("7. estimate iv on exposure beta and proportion of exposure variance explained") }
  IV_on_Exposure = est_beta_varexp( wdata = mod_data,
                                             response_variable = exposure,
                                             instrument = instrument,
                                             covariates = covariates,
                                             weights = weights)
  ## add some name specificity to these stats
  names(IV_on_Exposure) = paste0("iv_on_exposure_", names(IV_on_Exposure))
  
  ######################
  ### IX. IV on outcome 
  ###      beta & variance explained
  ######################
  if(messages == TRUE){ message("8. estimate iv on outcome beta and proportion of exposure variance explained") }
  IV_on_Outcome = est_beta_varexp( wdata = mod_data,
                                               response_variable = outcome,
                                               instrument = instrument,
                                               covariates = covariates,
                                               weights = weights)
  ## add some name specificity to these stats
  names(IV_on_Outcome) = paste0("iv_on_outcome_", names(IV_on_Outcome))
  
  ####################################
  ## RESULTS OUT
  ####################################
  if(messages == TRUE){ message("8. compiling data to report") }
  names(outcome) = "outcome"
  names(exposure) = "exposure"
  names(instrument) = "instrument"
  names(rnt_outcome) = "outcome_rnt"
  
  out = c( outcome, exposure, instrument, 
           sample_size,
           number_of_outcome_outliers, number_of_exposure_outliers, 
           rnt_outcome,
           W_outcome, W_exposure, 
           obs_stats , 
           IV_on_Exposure , 
           IV_on_Outcome ,
           mr_stats )
  
  
  if(messages == TRUE){ message("9. returning results to user") }
  return(out)
  
} ## end of function

obsmr_2_longformat = function(data, pop_id = "pop", alpha = 0.05, sandwich_estimates = FALSE){
  if(sandwich_estimates == FALSE){
    ## Observational
    a = data[, c("Obs_outcome","Obs_n","Obs_beta","Obs_se","Obs_P")]
    a$analysis = paste0("obs_", pop_id  )
    colnames(a) = c("outcome","n","beta","se","pval", "analysis")
    ## MR
    b = data[, c("MR_outcome","MR_n","MR_beta","MR_se","MR_P")]
    b$analysis = paste0("tsls_", pop_id  )
    colnames(b) = c("outcome","n","beta","se","pval", "analysis")
    
  }
  
  if(sandwich_estimates == TRUE){
    ## Observational
    a = data[, c("Obs_outcome","Obs_n","Obs_sw_beta","Obs_sw_se","Obs_sw_P")]
    a$analysis = paste0("obs_", pop_id  )
    colnames(a) = c("outcome","n","beta","se","pval", "analysis")
    ## MR
    b = data[, c("MR_outcome","MR_n","MR_sw_beta","MR_sw_se","MR_sw_P")]
    b$analysis = paste0("tsls_", pop_id  )
    colnames(b) = c("outcome","n","beta","se","pval", "analysis")
    
  }
  
  
  ## Define plot data
  out = rbind(a,b)
  ## insure it is a data frame
  out = as.data.frame(out)
  for(i in 2:5){out[,i] = as.numeric(out[,i])}
  
  ## add CIs
  out$lowerCI = out$beta - (1.96*out$se)
  out$upperCI = out$beta + (1.96*out$se)
  
  ## add a significance flag
  out$sig = "not_sig"
  w = which(out$pval <= alpha)
  out$sig[w] = "sig"
  
  ## set analysis levels
  out$analysis = as.factor(out$analysis)
  out$analysis = factor(out$analysis, levels = c( paste0("tsls_", pop_id ),
                                                  paste0("obs_", pop_id ) ) )
  
  out = out[, c("analysis","outcome","sig","n","beta","se","lowerCI","upperCI", "pval")]
  
  ## return
  return(out)
}


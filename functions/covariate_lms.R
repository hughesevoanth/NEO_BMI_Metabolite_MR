covariate_lms = function(datain, dependent, covariates, weights = NA ){
  data_out = sapply(covariates, function(cf){
    
    ## model formula
    form = as.formula( paste0( dependent, " ~", cf) )  
    
    ## linear model fitting
    if(is.na(weights[1])){
      fit = lm( form, data = datain )
    } else {
      fit = lm( form, data = datain, weights = weights  )
    }
    
    ## summary statistics
    s = summary(fit)
    n = length(residuals(s))
    a = anova(fit)
    etasq = a[1,2] / sum(a[,2])
    
    ## data to return
    out = c(dependent, cf, n, etasq, a[1,5])
    names(out) = c("dependent","covariate","sample_size","etasq","pvalue")
    return(out)
  })
  ###
  data_out = as.data.frame( t(data_out) )
  for(i in 3:5){data_out[,i] = as.numeric(data_out[,i])}
  ###
  return(data_out)

}

## Identify any IQR outliers and turn them into NAs. 
## outlier_values is the output from est_iqr(). 
## values is the vector of values that went into est_iqr().
na_iqr_outliers = function(values, outlier_values){
  neg = which( values < outlier_values[1] )
  pos = which( values > outlier_values[2] )
  outliers = c(neg, pos)
  if(length(outliers)>0){
    values[outliers] = NA
  }  
  return(values)
}

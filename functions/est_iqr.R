## Estimate cut off values for a vector (distribution) as the X ('distance') inter-quartile unit distances from the median.
est_iqr = function(x, distance = 5){
  m = median(x, na.rm = TRUE)
  q = quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr = c( m - (distance * q[1]),  
           m + (distance * q[2]) 
           )
  names(iqr) = c("lowerIQR","upperIQR")
  return( iqr )
}

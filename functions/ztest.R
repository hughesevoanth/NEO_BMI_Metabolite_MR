ztest = function(beta1, se1, beta2, se2){
  SE = as.numeric( sqrt( se1^2 + se2^2  ) )
  z = as.numeric( (beta2-beta1) / SE )
  P1 = pnorm(z, lower.tail = TRUE)  
  P2 = pnorm(z, lower.tail = FALSE)  
  P = min(P1,P2)
  out = c(z,P)
  return(out)
} 
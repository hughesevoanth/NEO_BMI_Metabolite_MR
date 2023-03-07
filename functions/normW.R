normW = function( y ){
  ## verify that x is a vector.
  if(is.vector(y) == FALSE){ stop("id_outliers() parameter y is expected to be a vector of data.") }

  ## verify that there is some variability
  if (length(unique(y)) == 1)
    stop("trait is monomorphic")
  if (length(unique(y)) == 2)
    stop("trait is binary")

  ## if there are more than 5000 observations
  ## sample down to n=5000
  if(length(y) <= 5000){
    W = shapiro.test( y )$stat; names(W) = "W"
  } else {
    W = shapiro.test( sample(y,5000) )$stat; names(W) = "W"
  }

  return(W)
}

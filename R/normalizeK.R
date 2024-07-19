normalizeK <- function(K){

  k = as.matrix(1/sqrt(diag(K)))
  nK = K * (k %*% t(k))
  
  return(nK)
  
}
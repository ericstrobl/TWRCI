remove_multicollinear <- function(X){
  
  tmp <- cor(X)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  
  keep = which(!apply(tmp, 2, function(x) any(abs(x) > 0.99, na.rm = TRUE)))
  
  # qr.X <- qr(X, tol=1e-9, LAPACK = FALSE)
  # rnkX <- qr.X$rank  ## 4 (number of non-collinear columns)
  # keep <- qr.X$pivot[seq_len(rnkX)]
  # ## 1 2 4 5 
  
  return(keep)
  
}
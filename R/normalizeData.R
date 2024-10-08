normalizeData <- function(X){
  
  X = as.matrix(X)
  for (i in seq_len(ncol(X))){
    if (sd(X[,i]) == 0){
      X[,i] = X[,i]-mean(X[,i])
    } else{
      X[,i] = (X[,i]-mean(X[,i]))/sd(X[,i]) # mean 1, sd 1
    }
  }
  X = as.matrix(X)
  
  return(X)
}

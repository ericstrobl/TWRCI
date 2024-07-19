ma <- function(x, n = 10){
  ma1 = filter(x, rep(1 / n, n), sides = 2)
  
  return(ma1)
  }



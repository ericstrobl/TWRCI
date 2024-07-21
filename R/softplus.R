softplus <- function(x,off=1){
  
  return( pmax(log(off),x)+log(1+exp(-abs(x-log(off)))) )
  
}
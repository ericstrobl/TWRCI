extract_causal_order <- function (G){
  require(MASS)
  p = ncol(G)
  iG = (ginv(diag(p)-G*0.5) > 1E-8)+0
  visit = sample(1:p,p)
  iG = iG[visit,visit]
  
  ord = which(rowSums(iG) == min(rowSums(iG)))[1]
  rem = setdiff(1:p,ord)
  while (length(rem)>1){
    ordn = which(rowSums(iG[-ord,-ord]) == min(rowSums(iG[-ord,-ord])))
    ord = c(rem[ordn[1]],ord)
    rem = setdiff(1:p,ord)
  }
  ord = c(rem,ord)
  
  return(visit[ord])
  
}
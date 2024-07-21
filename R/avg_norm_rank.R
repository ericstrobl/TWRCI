avg_norm_rank <- function(alphas,G,SNPs,alphas_GT = NULL){
  
  if (is.null(alphas_GT)){
    alphas_GT = SNP_weights_GT(G,SNPs)
  }
  alphas_GT = abs(alphas_GT)
  alphas_GT_ord = rank(alphas_GT)
  
  ix = which(alphas_GT > 1E-10)
  
  norm = mean(alphas_GT_ord[ix])-mean(1:length(ix))
  
  alphas_ord = rank(abs(alphas))
  avg_norm = mean(alphas_ord[ix])-mean(1:length(ix))
  
  if (norm == 0){
    return(1)
  }
  
  avg_norm = avg_norm/norm
  
  return(avg_norm)
}


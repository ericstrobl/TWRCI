SNP_weights_GT_genes <- function(G,SNPs){
  
  p = ncol(G$weights)
  
  betas = matrix(0,length(unlist(SNPs)),p)
  
  i_weights_X = ginv(diag(ncol(G$weights))-G$weights)
  for (g in 1:p){
    weights_X = i_weights_X[,g]
    weights_X[abs(weights_X)<1E-10] = 0
    
    for (s in 1:length(SNPs)){
      betas[SNPs[[s]],g] = G$weights_SNPs[[s]] * weights_X[s]
    }
  }
  
  return(betas)
}
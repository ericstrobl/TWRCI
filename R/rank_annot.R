rank_annot <- function(SNPs, SNPt, SNP_data, gene_data, G, K=NULL){
  
  SNP_data = normalizeData(SNP_data)
  gene_data = normalizeData(gene_data)

  betas_GT = SNP_weights_GT_genes(G,SNPt)
  
  L = length(SNPs)
  avg_norm_ranks = c()
  for (s in 1:L){
    if (length(SNPs[[s]])==0){
      avg_norm_ranks = c(avg_norm_ranks, 0)
    } else{
      betas = rep(0,length(unlist(SNPt)))
      if (!is.null(K)){
        iK = which(K==s)
        betas[SNPs[[s]]] = CV_LRR( cbind(SNP_data[,SNPs[[s]]],gene_data[,K[seq_len(iK-1)]]), gene_data[,s])$betas[1:length(SNPs[[s]])]
      } else{
        betas[SNPs[[s]]] = CV_LRR( SNP_data[,SNPs[[s]]], gene_data[,s])$betas
      }
      avg_norm_ranks = c(avg_norm_ranks, avg_norm_rank(betas,G,SNPs$SNPs,betas_GT[,s]))
    }
  }
  
  
  return( mean(avg_norm_ranks) )
  
  
}
sample_data <- function(G, SNPs, SNP_data, p=10, nsamps = 200){
  
  p = ncol(G$graph)
  
  # ix = sample(1:nrow(SNP_data),nsamps,replace=TRUE)
  # SNP_data = SNP_data[ix,,drop=FALSE]
  
  genes = get_gene_data(SNP_data,SNPs,G)
  # genes = sample_DAG_again(SNP_data,SNPs,G)
  
  
  return(  list(X = genes$data[,-G$Y], Y = genes$data[,G$Y], SNP_data = SNP_data, batch = genes$batches)   )
}
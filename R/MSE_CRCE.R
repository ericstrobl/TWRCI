MSE_CRCE <- function(Gest, SNPs, 
                     SNP_d,gene_data,batches,Y,
                     SNP_d_te,gene_data_te,batches_te,
                     truth=NULL){
  
  SK = SNP_d %*% t(SNP_d)
  SK_te = SNP_d_te %*% t(SNP_d_te)
  
  n = nrow(gene_data)
  
  cors = c()
  m0 = matrix(0,n,ncol(Gest)-1)
  for (k in 1:(ncol(Gest)-1)){
    
    # print(k)
    
    pa = which(Gest[,k])
    
    SK_i =  SNP_d[,SNPs[[k]],drop=FALSE] %*% t(SNP_d[,SNPs[[k]],drop=FALSE])
    SK_i = normalizeK(SK - SK_i)
    
    SK_i_te =  SNP_d_te[,SNPs[[k]],drop=FALSE] %*% t(SNP_d_te[,SNPs[[k]],drop=FALSE])
    SK_i_te = normalizeK(SK_te - SK_i_te)
    
    m1 = CV_KRR( cbind(gene_data[,pa],batches),SK_i, Y, cbind(gene_data_te[,pa],batches_te), SK_i_te)
    m2 = CV_KRR( cbind(gene_data[,c(pa,k)],batches),SK_i, Y, cbind(gene_data_te[,c(pa,k)],batches_te), SK_i_te)
    m0[,k] = m2-m1
    
  }
  
  if (is.null(truth)){
    return( m0 )
  } else{
    cors= c()
    for (i in 1:ncol(m0)){
      cors = c(cors, cor(m0[,i],truth[,i]))
      # cors = c(cors, sqrt(mean( (m0[,i] - truth[,i])^2 )) )
    }
    return( mean(cors) )
  }
  
}

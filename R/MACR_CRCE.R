MACR_CRCE <- function(Gest, SNPs,SNP_d,gene_data,batches,Y,
                                         SNP_d_te,gene_data_te,batches_te,Y_te){
  
  SK = rbind(SNP_d,SNP_d_te) %*% t(rbind(SNP_d,SNP_d_te))
  
  n = nrow(gene_data)
  
  cors = c()
  for (k in 1:(ncol(Gest)-1)){
    
    pa = which(Gest[,k])
    
    if ( (length(SNPs[[k]]) == 0) & (length(pa)==0) ){
      cors = c(cors, 1) # maximum penalty, because this should not be detected
    } else{
      
      if (length(SNPs[[k]])==0){next}
      
      SK_i =  rbind(SNP_d[,SNPs[[k]],drop=FALSE], SNP_d_te[,SNPs[[k]],drop=FALSE]) %*% t(rbind(SNP_d[,SNPs[[k]],drop=FALSE], SNP_d_te[,SNPs[[k]],drop=FALSE]))
      SK_i = normalizeK(SK - SK_i)
      
      m1 = CV_KRR( cbind(gene_data[,c(k,pa)],batches),SK_i[1:n,1:n], Y, 
                   cbind(gene_data_te[,c(k,pa)],batches_te), SK_i[(n+1):nrow(SK_i),1:n])
      
      res = Y_te - m1
      
      cors = c(cors, mean(abs(cor(res,SNP_d_te[,SNPs[[k]]]))) )
    }
    
  }
  
  print(cors)
  
  # plot(cors)
  
  return(mean(cors))
  
}

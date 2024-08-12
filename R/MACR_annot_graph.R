MACR_annot_graph <- function(G,SNPs,SNP_data,gene_data,batches,
                                         SNP_data_te,gene_data_te,batches_te,tar=1:ncol(G)){
  require(MASS)
  require(Rfast)
  require(earth)
  
  # gene_data = normalizeData(gene_data)
  
  n = nrow(gene_data)
  cors = c()
  d = length(SNPs)
  SNPt = 1:d
  for (k in tar){
    pa = which(G[,k])
    
    if (k>1){
      datat = cbind(SNP_data[,SNPs[[k]]],gene_data[,pa,drop=FALSE],batches)
      Xtet = cbind(SNP_data_te[,SNPs[[k]]],gene_data_te[,pa,drop=FALSE],batches_te)
    } else{
      datat = cbind(SNP_data[,SNPs[[k]]],batches)
      Xtet = cbind(SNP_data_te[,SNPs[[k]]],batches_te)
    }
    
    if (ncol(datat)==ncol(batches)){
      cors = c(cors, 1) # maximum penalty, because this should not be detected
    } else{
      # pre = CV_LRR2(datat,gene_data[,K[k],drop=FALSE],Xte=Xtet)$fitted.values
      # pre = predict(pre,Xtet)
      colnames(datat) = c()
      pre = earth(datat,gene_data[,k])
      # pre = predict(pre,Xtet)
      res = gene_data_te[,k]-predict(pre,Xtet)
      # if (length(unlist(SNPs))==0){
      #   cors = c(cors, 1)
      # } else{
      #   rem_SNPs = unlist(SNPs[K[(k+1):d]])
      #   cors = c(cors, mean(abs(cor(res,SNP_data_te[,rem_SNPs]))) )
      # }
      # pre = CV_LRR2(datat,gene_data[,K[k],drop=FALSE])$fitted.values
      # print(K[k])
      anc = setdiff(isAncAll(G,1:ncol(G),k),c(k,pa))
      if (length(anc)>0){
        geneM = mean(abs(cor(res,gene_data_te[,anc])))
        geneS = mean(abs(cor(res,SNP_data_te)))
        
        cors = c(cors, mean(c(geneM,geneS)) )
        # cors = c(cors, geneM )
      }
      # cors = c(cors, cor(pre,gene_data_te[,K[k],drop=FALSE]))
    }
    
  }
  
  # print(cors)
  
  # plot(cors)
  
  return(mean(cors))
  
}

get_genes_split <- function(SNP_data,ge_data,pvals,one.sided=TRUE){
  # partial correlation coefficient
  
  ge_data = normalizeData(ge_data)
  SNP_data = normalizeData(SNP_data)
  
  ## sample splitting hypothesis test
  groups = rep(1:2,length.out=nrow(ge_data))
  pM = c()
  
  for (g in 1:1){
    ig = which(groups==g)
    
    mod = CV_LRR2_pval(SNP_data[ig,,drop=FALSE],ge_data[ig,,drop=FALSE],
                  pvals,Xte=SNP_data[-ig,,drop=FALSE])
    # mod = CV_LRR2_pval(SNP_data[ig,mod$ix,drop=FALSE],ge_data[ig,,drop=FALSE],
    #                    pvals[mod$ix],Xte=SNP_data[-ig,mod$ix,drop=FALSE])
    
    ps = c()
    for (c in 1:ncol(ge_data)){
      if (one.sided){
        ps = c(ps, cor.test(ge_data[-ig,c],mod$fitted.values[,c],alternative="greater")$p.value)
      } else{
        ps = c(ps, cor.test(ge_data[-ig,c],mod$fitted.values[,c])$p.value)
      }
    }
    pM = cbind(pM,ps)
  }

  return(list(p=pM,ix = mod$ix))
  
  
}

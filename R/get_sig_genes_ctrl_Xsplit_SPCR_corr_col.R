get_sig_genes_ctrl_Xsplit_SPCR_corr_col <- function(SNP_data,ge_data,pvals,one.sided=TRUE){
  # partial correlation coefficient
  
  
  require(earth)
  
  ge_data = normalizeData(ge_data)
  SNP_data = normalizeData(SNP_data)
  
  ## sample splitting hypothesis test
  groups = rep(1:2,length.out=nrow(ge_data))
  
  ig = which(groups==1)
  
  mod1 = CV_LRR2_pval(SNP_data[ig,,drop=FALSE],ge_data[ig,,drop=FALSE],
                     pvals,Xte=SNP_data[-ig,,drop=FALSE])
  
  p = ncol(ge_data)
  
  ps = c()
  for (c in 1:ncol(ge_data)){
    if (one.sided){
      ps = c(ps, cor.test(ge_data[-ig,c],mod1$fitted.values[,c],alternative="greater")$p.value)
    } else{
      ps = c(ps, cor.test(ge_data[-ig,c],mod1$fitted.values[,c])$p.value)
    }
  }
  
  qs1 = qvalue(ps)$qvalues
  i1 = setdiff(which(qs1<0.1),p)
  print(i1)
  
  SNP_datat = SNP_data[, mod1$ix]
  for (s in 1:length(mod1$ix)){
    print(s)
    SNP_datat[,s] = earth(ge_data[,i1],SNP_data[, mod1$ix[s]])$residuals
  }
  
  ix = setdiff(1:length(ps),i1)
  ge_datat = ge_data[, -i1]
  for (t in 1:length(ix)){
    print(t)
    ge_datat[,t] = earth(ge_data[,i1],ge_data[,ix[t]])$residuals
  }
  
  mod2 = CV_LRR2_pval(SNP_datat[ig,,drop=FALSE],ge_datat[ig,,drop=FALSE],
                      pvals[mod1$ix],Xte=SNP_datat[-ig,,drop=FALSE],pt=1)
  
  ps = c()
  for (c in 1:ncol(ge_datat)){
    if (one.sided){
      ps = c(ps, cor.test(ge_datat[-ig,c],mod2$fitted.values[,c],alternative="greater")$p.value)
    } else{
      ps = c(ps, cor.test(ge_datat[-ig,c],mod2$fitted.values[,c])$p.value)
    }
  }
  
  qs2 = qvalue(ps)$qvalues
  i2 = setdiff(ix[which(qvalue(ps)$qvalues<0.1)],p)
  print(i2)

  
  return(list(i1 = i1, i2 = i2, qs1 = qs1, qs2 = qs2, ix = mod1$ix))
  
  
}

screen_anc_Y <- function (suffStat, Gest, SNPs=NULL, indepTest=earth_wrap, alpha=0.05, root=FALSE){
  
  p = ncol(Gest)
  # n = nrow(suffStat$data)
  # 
  # r = cor(suffStat$data[,p],suffStat$SNP_data)
  # t = r*sqrt((n-2)/(1-r^2))
  # iY = which((1-pt(abs(t),n-2))*2 < alpha)
  
  # vecs = vector("list",p-1)
  # for (i in 1:(p-1)){
  #   vecs[[i]] = hd.eigen(SNP_data_f[,aa$SNPs[[10]]],k=10,vectors=TRUE)$vectors
  # }
  # suffStat$batches = c()
  
  anc = c()
  for (x in 1:(p-1)){
    if (!is.null(SNPs)){
      if (root){
        pval <- indepTest(x, p, which(Gest[,x]), unlist(SNPs[-x]), suffStat)
        #print(pval)
      } else{
        pval <- indepTest(x, p, which(Gest[,x]), SNPs[[x]], suffStat)
        #print(pval)
      }
     
    } else{
      pval <- indepTest(x, p, which(Gest[,x]), NULL, suffStat)
    }
    
    if (pval < alpha){
      anc = c(anc, x)
    }
    
  }
  Gest[,p]=FALSE
  Gest[anc,p]=TRUE
  
  return(Gest)
}

susie_ss <- function(SNPs,betaY,betaY_var,cis,CL){

  # cis = get_cis_SNPs_window(locs, leads5, window=100000)
  
  require(coloc)
  
  # Y = normalizeData(Y)
  SNPs = normalizeData(SNPs)
  
  causalSNPs = c()
  for (i in 1:(length(cis)-1)){
    
    if (length(cis[[i]])==0){
      next
    } else if (length(cis[[i]])==1){
      causalSNPs <-   c(causalSNPs, cis[[i]])
      next
    }
    
    DY = list()
    DY$beta = betaY[cis[[i]]]
    DY$varbeta = betaY_var[cis[[i]]]
    DY$snp = as.character(cis[[i]])
    DY$position =  cis[[i]]
    DY$type = "quant"
    DY$sdY = 1

    DY$LD = cor(SNPs[,cis[[i]],drop=FALSE])
    colnames(DY$LD) = DY$snp; rownames(DY$LD) = DY$snp
    DY$N = nrow(SNPs)

    causalSNPs <-   c(causalSNPs, cis[[i]][runsusie(DY)$pip>0.8])
  }
  
  alpha = rep(0,ncol(SNPs))
  n = nrow(SNPs)
  betaY = as.matrix(betaY)
  if (length(causalSNPs)>0){
    alpha[causalSNPs] = ginv(cor(SNPs[,causalSNPs,drop=FALSE])*n) %*% (betaY[causalSNPs,drop=FALSE]*n)
  }

  alpha_max = rep(0,length(CL))
  for (c in 1:length(CL)){
    alpha_max[c] = max(abs(alpha[CL[[c]]]))
  }
  alpha_max[is.infinite(alpha_max)]=0
  
  return(list(alpha=alpha,alpha_max=alpha_max))
  
  
  
  
}
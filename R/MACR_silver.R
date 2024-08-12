MACR_silver <- function(Y,G,exp_d,SNP_d,SNPs,nu_d,desc,exp_d_all){
  
  anc = isAncAll(G,1:(ncol(G)-1),ncol(G))
  
  desc=setdiff(desc,anc)
  
  resids1 = earth( cbind(exp_d[,anc],SNP_d[,SNPs[[ncol(G)]]],nu_d),Y)$residuals
  resids2 = earth( cbind(exp_d[,anc],SNP_d[,SNPs[[ncol(G)]]],nu_d),exp_d_all[,desc])$residuals

  return( mean(abs(cor(resids1,resids2))) )
}

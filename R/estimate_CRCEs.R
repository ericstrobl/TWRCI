estimate_CRCEs <- function(X,batches,Gest,SNPs,SNP_d,Y,anc,Xte=NULL){
  
  if (!is.null(Xte)){
    Xte = X
    SNP_d_te = SNP_d
  }

  CRCE = c()
  CRCE_mean = c()
  p = ncol(X)+1
  
  SNP_d=normalizeData(SNP_d)
  
  SK = SNP_d %*% t(SNP_d)
  
  for (i in anc){
    # anc = K[1:which(K==i)]
    anc = c(which(Gest[,i]),i)

    if (is.null(Xte)){
      SNP_d_i = SNP_d[,SNPs[[i]]]
      SK_i =  SNP_d_i %*% t(SNP_d_i)
      SK_i = normalizeK(SK - SK_i)

      m1 = CV_KRR( cbind(X[,anc],batches),SK_i, Y )

      m2 = CV_KRR( cbind(X[,anc[-length(anc)]],batches), SK_i, Y)

    } else{
      # m1 = CV_KRR( cbind(X[,anc],SNP_d[,unlist(SNPs[-i])]), Y,
      #              Xte = cbind(Xte[,anc], SNP_d_te[,unlist(SNPs[-i])]) )
      # 
      # m2 = CV_KRR( cbind(X[,anc[-length(anc)]], wS*SNP_d[,unlist(SNPs[-i])]), Y,
      #              Xte = cbind(Xte[,anc[-length(anc)]], wS*SNP_d_te[,unlist(SNPs[-i])]) )
    }
    # 
    CRCE_mean = c(CRCE_mean,  mean( m1[Y>1] - m2[Y>1]) - mean( m1[Y<=1] - m2[Y<=1]))
    CRCE = cbind(CRCE,  m1 - m2)
    # print(CRCE_mean)
  }
  
  return( list(CRCE_mean=CRCE_mean, CRCE = CRCE))
  
  
}

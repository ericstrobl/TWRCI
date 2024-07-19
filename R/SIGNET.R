SIGNET <- function(X,batches,SNPd,cis){
  
  X = earth(batches,X)$residuals
  
  require(glmnet)
  d = ncol(X)
  # cis = get_cis_SNPs(locs,ncol(X),ncol(SNPd))
  # cis = locs
  
  Xe = X
  for (i in 1:d){
    print(i)
    Xe[,i] = CV_LRR(SNPd,X[,i])$fitted.values
  }
  
  ix = 1:d
  G = matrix(0,d,d)
  for (i in 1:d){
    # adaptive lasso
    feat = cbind(Xe[,-i],SNPd[,cis[[i]]])
      
    betas = CV_LRR(feat, X[,i])$betas
    
    feat = sweep(feat,2,1/betas,'*')
    
    cv_model <- cv.glmnet(feat, X[,i], alpha = 1)
    best_model <- glmnet(feat, X[,i], alpha = 1, lambda = cv_model$lambda.min)
    
    pa = which( abs(coef(best_model)[2:d]/betas[1:(d-1)]) > 5E-2)
    pa = setdiff(1:d,i)[pa]
    G[pa,i] = 1
  }
  
  return(G)
  
}
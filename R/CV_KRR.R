CV_KRR <- function(X,SK,Y,Xte=NULL,SKte=NULL){
  require(Rfast)
  Y=as.matrix(Y)
  
  d=ncol(Y)
  mY = colMeans(Y)
  sY = apply(Y,2,sd)
  for (i in 1:d){
    Y[,i] = (Y[,i,drop=FALSE] - mY[i])/sY[i]
  }
  
  n = nrow(Y);
  
  # if (is.null(Xte)){
  #   X=normalizeData(X)
  # } else{
  #   Xte=as.matrix(Xte)
  #   Xp=normalizeData(rbind(X,Xte))
  #   X = Xp[1:n,,drop=FALSE]
  #   Xte = Xp[(n+1):nrow(Xp),,drop=FALSE]
  # }
  
  # shuf = sample(1:n,n,replace=FALSE)
  # X[1:n,]=X[shuf,]
  # Y[1:n,] = Y[shuf,]
  if (length(X)>0){
    X=as.matrix(X)
    dotX = as.matrix(proxy:::dist(rbind(X,Xte),y=X))
    med_dist = median(dotX);
    if (med_dist == 0){ med_dist = 1;}
    # med_dist = 1
    dotX = dotX^2;
    sigmas = med_dist*seq(0.5,4,0.5);
  } else{
    sigmas = 1
  }
  
  lambdas = c(100,10,1,1E-1,1E-2)
  
  ## estimate conditional expectations
  err = Inf; #maximum error allowed
  
  for (s in seq_len(length(sigmas))){
    if (length(X)>0){
      K = RBF_kernel(dotX[1:n,],sigmas[s]) + SK
    } else{
      K = SK
    }
    # Ke = eigen(K,symmetric=TRUE)
    for (l in seq_len(length(lambdas))){
      # K_inv = Ke$vector %*% diag(1/(Ke$values+n*lambdas[l])) %*% t(Ke$vector)
      K_inv = spdinv( K  + n*diag(n)*lambdas[l]);
      H = K_inv %*% K;
      te_y = t(t(Y) %*% H)
      res = LOO_KRR(Y, Y - te_y, diag(H))
      if (res$err<err){
        err = res$err
        pre = res$pre
        sigma_f = sigmas[s]
        alpha = K_inv %*% Y
      }
    }
  }
  
  if (!is.null(Xte)){
    Kte = RBF_kernel(dotX[(n+1):nrow(dotX),],sigma_f)
    pre = (Kte + SKte) %*% alpha
  }
  
  pre = pre*sY[i] + mY[i]
  
  return(pre)
  
}
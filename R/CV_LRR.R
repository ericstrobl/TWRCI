CV_LRR <- function(X,Y,maxY=NULL,minY=NULL,Ci=NULL,lambda=NULL,Xte=NULL){
  require(Rfast)
  
  Y = as.matrix(Y)
  dY = ncol(Y)
  n = nrow(Y)
  
  if ( length(X)==0 ){
    return(list(betas = 0, lambda = 0, fitted.values = t(replicate(n, colMeans(Y))) ))
  }
  
  X = as.matrix(X)
  d = ncol(X)
  
  
  dotX = t(X) %*% X
  dotXY = t(X) %*% Y
  
  # print(length(X))
  # print(nrow(X))
  # print(ncol(X))
  if (is.null(lambda)){
    lambdas = c(100,10,1,1E-1,1E-2,1E-3);
    # lambdas = c(0.1,0.2,0.3,0.5,0.7,0.9)
  } else{
    lambdas = lambda
  }
  
  ## estimate conditional expectations
  err = rep(Inf,dY)
  lambdaf = rep(Inf,dY)
  beta = matrix(0,d,dY); 
  if (is.null(Xte)){
    pre = matrix(0,n,dY) 
  } else{
    pre = matrix(0,nrow(Xte),dY) 
  }
  for (l in seq_len(length(lambdas))){
    I = diag(d)
    if (!is.null(Ci)){
      I[Ci:ncol(X),Ci:ncol(X)]=0
    }
    U = chol( dotX  + n*I*lambdas[l])
    dotX_inv = chol2inv(U);
    Z = forwardsolve(t(U),t(X))
    dH = colSums(Z^2)
    # H = X %*% dotX_inv %*% t(X);
    te_y =  X%*% dotX_inv %*% dotXY
    res = LOO_KRR(Y, Y - te_y, dH)
    
    iy  = which(res$err<err)
    err[iy] = res$err[iy]
    lambdaf[iy] = lambdas[l]
    beta[,iy] = dotX_inv %*% dotXY[,iy]
    if (is.null(Xte)){
      pre[,iy] = te_y[,iy]
    } else{
      pre[,iy] = Xte %*% beta[,iy]
    }

  }
  
  return(list(betas = abs(beta), lambda = lambdaf, fitted.values = pre))
  
}
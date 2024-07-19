CV_LRR2_pval <- function(Xt,Y,pvals,lambda=NULL,Xte=NULL,pt=NULL){
  require(Rfast)
  
  Y = as.matrix(Y)
  dY = ncol(Y)
  n = nrow(Y)
  
  if ( length(X)==0 ){
    return(list(betas = 0, lambda = 0, fitted.values = t(replicate(n, colMeans(Y))) ))
  }
  
  if (is.null(lambda)){
    lambdas = 1 ###
  } else{
    lambdas = lambda
  }
  
  if (is.null(pt)){
    pt = quantile(pvals,seq(1,0.1,length.out=10))
  }
  dt = ncol(Xt)
  beta = matrix(0,dt,dY); 
  if (is.null(Xte)){
    pre = matrix(0,n,dY) 
  } else{
    pre = matrix(0,nrow(Xte),dY) 
  }
  
  err_overall = Inf
  lambdaf = rep(Inf,dY)
  for (pv in pt){
    err = rep(Inf,dY)
    ix = which(pvals<=pv)
    X = as.matrix(Xt[,ix])
    if (length(X)==0){
      next
    }
    d = ncol(X)
    
    if (n >= d){
      dotX = t(X) %*% X
      dotXY = t(X) %*% Y
    } else{
      dotX = X %*% t(X)
    }
    
    ## estimate conditional expectations
    
    for (l in seq_len(length(lambdas))){
      if (n>=d){
        I = diag(d)
        U = chol( dotX  + n*I*lambdas[l])
        dotX_inv = chol2inv(U);
        Z = forwardsolve(t(U),t(X))
        dH = colSums(Z^2)
        te_y =  X%*% dotX_inv %*% dotXY
        # print(max(dH))
      } else{
        I = diag(n)
        dotX_inv = spdinv( dotX  + n*I*lambdas[l])
        Z = dotX %*% dotX_inv
        dH = diag(Z)
        te_y =  Z %*% Y
        # print(max(dH))
      }
      res = LOO_KRR_cor(Y, Y - te_y, dH)
      
      iy  = which(res$err<err)
      err[iy] = res$err[iy]
      
      if (mean(err)<err_overall){
        err_overall=mean(err)
        
        pf = pv
        ixf = ix
        
        lambdaf[iy] = lambdas[l]
        if (n>=d){
          beta[ixf,iy] = dotX_inv %*% dotXY[,iy,drop=FALSE]
        } else{
          beta[ixf,iy] = t(X) %*% dotX_inv %*% Y[,iy,drop=FALSE]
        }
        if (is.null(Xte)){
          pre[,iy] = te_y[,iy]
        } else{
          pre[,iy] = Xte[,ixf] %*% beta[ixf,iy]
        }
      }
      
    }
    
  }
  beta = beta[ixf,,drop=FALSE]
  
  return(list(betas = abs(beta), lambda = lambdaf, fitted.values = pre, thres = pf,
              ix = ixf, beta = beta))
  
}
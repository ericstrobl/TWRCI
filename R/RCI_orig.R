RCI_orig <- function(X,alpha=0.2){
  time = proc.time()
  K = DirectLiNGAM_fast_orig(X,Y) # Local Plus
  return(K)
}

DirectLiNGAM_fast_orig <- function(X,alpha=0.2){
  X = as.matrix(X)
  
  K = c() #1 
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  U = 1:ncol(X) #2
  
  repeat{ #11
    
    # U = U[cortest_fast(X[,U,drop=FALSE],Y,alpha=alpha)] ###
    # if (length(U)==0){ ###
    #   break ###
    # } ###

    root = FindRoot_fast_orig(X,U,Cov) #6
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    if (length(U)==0){ ###
      break ###
    } ###
    X = UpdateData_orig(X,U,Cov,root) #9
    Cov = UpdateCovMat_orig(U,Cov,root) #10
  }
  
  return(K)
}

UpdateCovMat_orig <- function(U,Cov,root){
  
  for (j in setdiff(seq_len(ncol(Cov)),root)){
    Cov[U,j] = (Cov[U,j] - Cov[U,root] * Cov[root,j])/
      (sqrt(1-Cov[U,root]^2)*sqrt(1-Cov[j,root]^2)) 
  }
  
  return(Cov)
  
}

UpdateData_orig <- function(X,U,Cov,root){
  
  for (j in setdiff(U,root)){
    X[,j] = (X[,j] - Cov[j,root] / Cov[root,root] * X[,root,drop=FALSE])/
      sqrt(1-Cov[j,root]^2)
  }
  
  return(X)
  
}

FindRoot_fast_orig <- function(X,U,Cov){
  
  r = length(U) #4
  
  if (r==1){ #1
    return(U) #2
  }
  
  O = rep(1,r)
  S = rep(0,r)
  M = matrix(TRUE,r,r)
  
  repeat{
    I = (S == min(S))
    
    if ( sum(O[I] > r) > 0){
      break
    }
    
    for (i in which(I)){
      if (M[i,O[i]]){
        score = Compare2(X,U[i],U[O[i]],Cov)
        S[i] = S[i] + min(0,score)^2
        S[O[i]] = S[O[i]] + min(0,-score)^2
        
        M[i,O[i]]=FALSE
        M[O[i],i]=FALSE
      }
      
      O[i] = O[i]+1
    }
  }
  
  root = U[S==min(S)][1]
  
  return(root) #output
  
}

cortest_fast <- function(X,Y,suffStat=NULL,alpha=0.2){
  
  n = nrow(X)
  rs = cor(X,Y)
  t = sqrt(n-2)*rs/sqrt(1-rs^2); ps = 1-pf(t^2,1,n-2) # f test
  
  return( which(ps<alpha) )
  
}

Compare2 <- function(X,i,j,Cov){
  
  if (i==j){  ####
    return(0) ####
  }
  
  rij = X[,i] - Cov[i,j]/Cov[j,j] * X[,j] #9
  rji = X[,j] - Cov[j,i]/Cov[i,i] * X[,i] #10
  
  rsd = sqrt(1-Cov[i,j]^2)
  rij = rij / rsd #11
  rji = rji / rsd #12
  
  score = LRT(X[,i],X[,j],rij,rji)
  
  return(score) #13  ####
}

LRT <- function(xi,xj,rij,rji){
  
  return(Hu(xj)+Hu(rij)-Hu(xi)-Hu(rji))
}

Hu <- function(u){
  
  k1 = 79.047
  k2 = 7.4129
  beta = 0.37457
  
  H = 0.5*(1+log(2*pi))-k1*(mean(log(cosh(u)))-beta)^2 - k2*mean(u*exp(-(u^2)/2))^2
  
  return(H)
}
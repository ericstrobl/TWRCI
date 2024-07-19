GRCI_ANM <- function(X){
  require(RANN)
  require(gtools)
  require(earth)
  require(independence)
  
  p = ncol(X)
  G = matrix(TRUE,p,p)
  diag(G) = FALSE
  
  K = DirectHNM_fast_Y_ANM(X,G)
  
  return(rev(K))
}

DirectHNM_fast_Y_ANM <- function(X,G){
  X = as.matrix(X)
  # E = X
  
  K = c()
  S = rep(-Inf,ncol(X))
  U = 1:ncol(X)
  update = U
  
  repeat{ 
    s_out = FindSink_ANM(X,U,S,G,update)
    sink = s_out$sink
    S = s_out$S
    
    K = c(K,sink)
    U = U[-which(U==sink)]
    
    if (length(U)==0){ ###
      break ###
    } ###
    
    if (sum(G[sink,])){
      # X[,sink] = PartialOut(X[,G[sink,]],X[,sink],w[G[sink,],,drop=FALSE]) # partial out neighbors from sink
      # X[,sink] = nb_Pearson_ANM_K08(X[,G[sink,],drop=FALSE],X[,sink],C,Xp)
      X[,sink] = earth(X[,G[sink,],drop=FALSE],X[,sink])$residuals
    }
    update = intersect(U,which(G[sink,]))
    G[sink,] = FALSE; G[,sink]=FALSE # remove node
    
  }
  
  return(K) #output
}


FindSink_ANM <- function(X,U,S,G,update,w=NULL){
  r = length(update)
  
  ## BASELINE SCORES
  for (i in seq_len(r)){
    V = which(G[update[i],])
    if (length(V)==0){ # if no neighbors, then break because score is -Inf
      next
    }
    
    S[update[i]] = CompareG_ANM(X,update[i],V,w) # baseline
    
  }
  
  sink = U[S[U]==min(S[U])][1] 
  
  return(list(sink=sink,S=S)) #output
  
}

CompareG_ANM <- function(X,i,js,w=NULL){
  
  if (identical(i,js)){  ####
    return(0) ####
  }
  
  rXY = earth(X[,js,drop=FALSE],X[,i])$residuals
  
  score = c()
  for (k in 1:length(js)){
    score = c(score, hoeffding.D.test(X[,js[k]]+rnorm(length(X[,js[k]]))*1E-6,
                                      c(rXY)+rnorm(length(rXY))*1E-6,precision=1)$scale)
  }
  score = max(score)
  
  return(score) #13  ####
}

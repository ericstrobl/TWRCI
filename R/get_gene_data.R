get_gene_data <- function(SNP_d,SNPs,DAG){
  
  # SNP_d = normalizeData(SNP_d)
  require(MASS)
  
  G = DAG$graph
  p = nrow(G)
  
  nsamps = nrow(SNP_d)
  
  dataT = matrix(0,nsamps,p)
  dataE = dataT
  off = DAG$off
  
  L = length(DAG$ord)
  for (i in 1:L){
    c = DAG$ord[i]
    h = which(G[DAG$ord[1:(i-1)],c]==1)
    h = DAG$ord[1:(i-1)][h]
    
    if (length(h)>0){
      A = dataT[,h,drop=FALSE]%*%as.matrix(DAG$weights[h,c,drop=FALSE]) # do not need to include offset because the gamma error terms introduce the offset
    } else{
      A = rep(0,nsamps)
    }
    
    dataE[,c] = 0.2*rnorm(nsamps)
    dataT[,c]= A + dataE[,c] + SNP_d[,SNPs[[c]]]%*%DAG$weights_SNPs[[c]]
    
  }

  # data = dataT
  # batches = rep(1,nsamps)

  # dataT = cbind(softplus(dataT[,-DAG$Y]),dataT[,DAG$Y])
  dataT = softplus(dataT)

  # batch modification
  batches = sample(DAG$n_batch,nsamps,replace=TRUE) # each sample is from the same batch
  datap = dataT
  for (b in 1:DAG$n_batch){
    ib = which(batches == b)
    datap[ib,-DAG$Y] = t(t(dataT[ib,-DAG$Y])* DAG$batch[b,-DAG$Y])  #####
  }
  # # introduce Poisson measurement error from multinomial
  total = rowSums(datap[,-DAG$Y])###

  datap[,-DAG$Y] = datap[,-DAG$Y]/total # compute proportions ###

  N = rep(0,nsamps)
  data = dataT
  for (n in 1:nsamps){ # sample from multinomial with proportions
    N[n] = rpois(1,total[n])
    data[n,-DAG$Y] = rmultinom(1,N[n],prob=datap[n,-DAG$Y])
  } # measurement error marginal distributions then follow a Poisson


  return(list(data=data,dataT=dataT,dataE=dataE,batches=batches))

  # return(  list(data=data,dataT = dataT,batches=batches,dataE=dataE)  )
}


generate_DAG <- function(p,SNPs,en=2,bulk_nbatch=1){
  
  require(Matrix)
  
  N = p*p-p;
  
  DAGs = list()
  DAGs$Y = 0
  DAGs$HK = 0
  
  ps = unlist(length(SNPs))
  graph = matrix(0,p,p)
  while( DAGs$Y != p  ){
    
    ### SINGLE CELL
    
    samplesB = rbinom(N/2,1, en/(p-1) ); # sample edges
    graph[upper.tri(graph, diag=FALSE)] <- samplesB; # put in edges in upper triangular
    
    ord = sample(1:p,p,replace=FALSE) # permute order
    DAGs$graph = graph[ord,ord]
    # ord2 = rep(0,p)
    # ord2[ord] = 1:p
    DAGs$ord = ord
    
    Ys = which( (rowSums(DAGs$graph)==0) & (colSums(DAGs$graph)>0) ) # no children, some observed parents
    
    if (p %in% Ys) { # ensure that DAGs$Y is the last variable
      DAGs$Y = p
    }
  }
  
  
  DAGs$graph = as.matrix(DAGs$graph)
  
  # partial order
  ord = rep(0,p)
  ord[DAGs$ord] = 1:p
  DAGs$ord = ord
  
  ix = which(DAGs$graph!=0)
  weights = matrix(0.75*runif(length(ix))+0.25)*sample(c(-1,1),length(ix),replace=TRUE)
  DAGs$weights = DAGs$graph
  DAGs$weights[ix] = weights
  
  weights_SNPs = vector("list",length(SNPs))
  for (s in 1:length(SNPs)){
    weights_SNPs[[s]] = matrix(0.1*runif(length(SNPs[[s]]))+0.05)*sample(c(-1,1),length(SNPs[[s]]),replace=TRUE)
  }
  DAGs$weights_SNPs = weights_SNPs
  
  DAGs$off = runif(p)*0.4 + 0.1 #

  DAGs$n_batch = bulk_nbatch
  DAGs$batch = runif(p*DAGs$n_batch)*9900 + 100 #10 to 200, volume times pi_{ij}
  DAGs$batch = matrix(DAGs$batch,DAGs$n_batch,p) #batches are rows, variables are columns
  
  return(DAGs)
}
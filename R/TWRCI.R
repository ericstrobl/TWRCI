TWRCI <- function(genes_d,SNPs_d,target_d,nuisance_d,batches){
  # genes_d = gene expression data
  # SNPs_d = variant data
  # target_d = phenotype data
  # nuisance_d = nuisance variable data
  # batches = batches
  
  require(earth)
  
  SNPs_d = normalizeData(SNPs_d)
  genes_d = normalizeData(genes_d)
  target_d = normalizeData(target_d)
  batches = normalizeData(batches)
  
  DL = DirectLiNGAM_SG_11_no_ss2(genes_d,SNPs_d,target_d,nuisance_d,batches)
  
  return(list(K = DL$K, SNPs = DL$SNPs))
}

DirectLiNGAM_SG_11_no_ss2 <- function(genes_d,SNPs_d,target_d,nuisance_d,batches){
  
  yi = ncol(genes_d)+1
  
  SNPs = vector("list",yi)
  
  #### for Y only
  LRR = CV_LRR2_all(cbind(SNPs_d,batches),genes_d)
  betas = LRR$betas
  lambda = LRR$lambda
  
  betas_i = CV_LRR2_all(cbind(SNPs_d,genes_d,batches),target_d,lambda=lambda)$betas # dont need nuisance_d here
  
  iSN = c()
  for (s in 1:ncol(SNPs_d)){ # find SNPs with largest coefficients
    if (betas_i[s] > max(betas[s,])){
      iSN = c(iSN,s)
    }
  }
  
  if (length(iSN)>0){
    SNPs[[yi]] = iSN
  } else{
    SNPs[yi] = list(NULL)
  }
  
  K = yi
  U = 1:(yi-1)
  
  repeat{#11
    
    FS = FindSink_GS_112(U,genes_d,SNPs_d,nuisance_d,batches,SNPs) #6
    sink = FS$sink
    K = c(sink,K) #7
    U = U[-which(U==sink)] #8
    
    if (length(FS$SNP)>0){
      SNPs[[sink]] = FS$SNP
    } else{
      SNPs[sink] = list(NULL)
    }
    
    if (length(U)==0){ ###
      break ###
    } ###
    
  }
  
  
  return(list(K=K, SNPs=SNPs) ) #output
}

FindSink_GS_112 <- function(U,genes_d,SNPs_d,nuisance_d,batches,SNPs){
  
  r = length(U) #4
  id = setdiff(1:ncol(SNPs_d),unlist(SNPs))
  if (r==1){
    return(list(sink=U, SNP = id))
  }
  
  LRR = CV_LRR2_all(cbind(SNPs_d[,id,drop=FALSE],batches),genes_d[,U,drop=FALSE])
  betas = LRR$betas
  lambda = LRR$lambda
  
  
  ## BASELINE SCORES
  stat = Inf
  for (i in seq_len(r)){
    PS = Compare_ps_sink_112(U[i],U[-i],betas[,-i,drop=FALSE],lambda,genes_d,SNPs_d,nuisance_d,batches,SNPs)
    
    if (PS$stat <= stat){
      stat = PS$stat
      sink = U[i]
      SNP = PS$SNP
    }
    
  }
  
  return(list(sink=sink, SNP = SNP)) #output
  
}

Compare_ps_sink_112 <- function(i,not_is,betas,lambda,genes_d,SNPs_d,nuisance_d,batches,SNPs){
  n = nrow(genes_d)
  id = setdiff(1:ncol(SNPs_d),unlist(SNPs))
  if (length(id)==0){
    return(list(stat = Inf, SNP = NULL))
  }
  
  # betas_i = earth_beta(cbind(SNPs_d[,id,drop=FALSE],genes_d[,not_is]),genes_d[,i])
  betas_i = CV_LRR2_all(cbind(SNPs_d[,id,drop=FALSE],genes_d[,not_is],nuisance_d,batches),genes_d[,i],lambda=lambda)$betas
  
  iSN = c()
  for (s in 1:length(id)){ # find SNPs with largest coefficients
    if (betas_i[s] > max(betas[s,])){
      iSN = c(iSN,id[s])
    }
  }
  
  if ( (length(iSN)==0) | (length(setdiff(id,iSN))==0) ){
    return(list(stat = Inf, SNP = iSN))
  }
  
  # print(length(colnames(genes_d[,not_is])))
  # print(length(unique(colnames(genes_d[,not_is]))))
  
  # print(length(colnames(SNPs_d[,iSN])))
  # print(length(unique(colnames(SNPs_d[,iSN]))))
  
  resid = earth(cbind(genes_d[,not_is],SNPs_d[,iSN],nuisance_d,batches),genes_d[,i])$residuals
  
  rs = cor(resid,SNPs_d[,setdiff(id,iSN)]) # correlate to other SNPs
  # rs = cor(resid,SNPs_d)
  
  return(list(stat = mean(abs(rs)), SNP = iSN))
  
}

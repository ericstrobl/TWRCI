get_coloc_ss <- function(X,SNPs,batches,betaY,betaY_var,cis,method="abf"){
  
  require(coloc)
  
  X = normalizeData(X)
  # Y = normalizeData(Y)
  SNPs = normalizeData(SNPs)
  
  n = nrow(SNPs); d = ncol(SNPs)
  PP4 = matrix(0,d,ncol(X))

  for (i in 1:ncol(X)){
    
    if (length(cis[[i]])==1){
      PP4[cis[[i]],i] = 1
      next
    } else if(length(cis[[i]])==0){
      next
    }
    
    # betaY = c()
    # betaY_var = c()
    # n = nrow(X)
    # for (s in cis[[i]]){
    #   mod = lm.fit(SNPs[,s,drop=FALSE],Y)
    #   betaY = c(betaY, mod$coefficients[1])
    #   betaY_var = c(betaY_var, var(mod$residuals)/(n-1))
    # }
    
    DY = list()
    DY$beta = betaY[cis[[i]]]
    DY$varbeta = betaY_var[cis[[i]]]
    DY$snp = as.character(cis[[i]])
    DY$position =  cis[[i]]
    DY$type = "quant"
    DY$sdY = 1
    
    n = nrow(X)
    prod = SNPs[,cis[[i]],drop=FALSE]*c(earth(batches,X[,i])$residuals)
    betaX = colMeans(prod)
    betaX_var = apply(prod,2,var)/n
    
    DX = list()
    DX$beta = betaX
    DX$varbeta = betaX_var
    DX$snp = as.character(cis[[i]])
    DX$position = cis[[i]]
    DX$type = "quant"
    DX$sdY = 1
    
    if (method == "susie"){
      DY$LD = DX$LD = cor(SNPs[,cis[[i]],drop=FALSE])
      colnames(DY$LD) = DY$snp; rownames(DY$LD) = DY$snp
      colnames(DX$LD) = DX$snp; rownames(DX$LD) = DX$snp
      DY$N = DX$N = nrow(SNPs)
      pos <- coloc.susie(dataset1=DY, dataset2=DX)$results[,2]
      if (length(pos)==0){
        pos = 0
      } else{
        pos = as.matrix(pos)
      }
    } else{
      pos <- coloc.abf(dataset1=DY, dataset2=DX)$results[,12]
    }
    PP4[cis[[i]],i] = pos
    
  }
  
  SNPs = vector("list",ncol(X)+1)
  
  for (i in 1:ncol(X)){
    SNPs[[i]] = which(PP4[,i] > apply(PP4[,-i],1,max))
  }
  
  
  return(SNPs)
  
  
  
  
}
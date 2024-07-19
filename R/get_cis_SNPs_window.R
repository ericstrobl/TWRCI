get_cis_SNPs_window <- function(locs, leads5, window=500000, overlap=FALSE){
  
  require(FNN)
  
  leads_locs = leads5[unlist(locs),]
  
  SNPs = list()
  
  SNP_c = c(); is = c(); distances = c()
  for (i in 1:nrow(leads_locs)){
    SNPs[[i]] = which( (leads5$chr==leads_locs$chr[i]) &  # includes at least itself
                     (leads5$position <= (leads_locs$position[i] + window) ) &
                     (leads5$position >= (leads_locs$position[i] - window) ) )
    
    if (overlap == FALSE){
      SNP_c = c (SNP_c,SNPs[[i]] )
      is = c(is, rep(i,length(SNPs[[i]])))
      distances = c(distances, abs(leads5$position[SNPs[[i]]] - leads_locs$position[i])) 
    }
    
  }
  
  if (overlap == FALSE){
    or = order(distances)
    ix = which(!duplicated(SNP_c[or]))
    SNP_c = SNP_c[or][ix]
    is = is[or][ix]
    
    SNPs = list()
    for (i in 1:nrow(leads_locs) ){
      SNPs[[i]] = sort(SNP_c[is == i])
    }
  }
  
  return(SNPs)
  
}

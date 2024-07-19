generate_SNP_list <- function(p,q,leads){
  # p = number of genes
  # q = number of SNPs
  
  mult = q/p
  
  SNP_count = 1+rmultinom(n = 1, size = round(mult*p)-p, prob = rep(1/p, p))
  SNPs = list()
  
  
  cand = 1:q
  locs = as.list(sample(cand,p,replace=FALSE)) # random locations to assess algorithms above and beyond prior knowledge
  for (i in sample(1:p,p)){
    
    # prob = locs[[i]]/(q+1)
    # probs = pmax(dbinom(0:(q-1),q,prob),1/q)
    # probs = probs/sum(probs)
    
    # print(i)
    intra = which(leads$chr[cand]==leads$chr[locs[[i]]])
     # farther from intra-chromosone
    if (length(intra)>0){
      outer = setdiff(1:length(cand),intra)
      
      probs = rep(0,length(cand))
      probs[intra] = 1/(1+order(abs(leads$position[cand][intra] - leads$position[locs[[i]]])))
      probs[outer] = 1/(1+max(order(abs(leads$position[cand][intra] - leads$position[locs[[i]]])))+1)
      probs = probs/sum(probs)
    } else{
      probs = rep(1/length(cand),length(cand))
    }
    
    # probs = rep(1/q,q)
    
    if (length(cand)>1){
      SNPs[[i]] = sample(cand,SNP_count[i],pmax(probs,1E-200),replace=FALSE)
    } else{
      SNPs[[i]] = cand
    }
    
    cand = setdiff(cand,SNPs[[i]])
  }
  
  return(list(SNPs=SNPs,locs=locs))
  
}
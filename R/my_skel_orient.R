my_skel_orient <- function(suffStat, K, SNPs=NULL, alpha=0.01){
  
  p = length(K)
  
  G_est1 = my_skeleton_p(suffStat, p, SNPs=SNPs, alpha=alpha)$G
  for (k in seq_len(p)){
    G_est1[K[k],K[seq_len(k-1)]]=FALSE
  }
  # 
  # return(G_est1)
  G_est = G_est1
  G_est[,p]=TRUE
  G_est2 = my_skeleton_p(suffStat, p, SNPs=SNPs, Gest = G_est, alpha=min(0.05,alpha*5))$G # check graph consistency with binary target
  
  return(list(G_est2 = G_est2,G_est1 = G_est1))
  
}
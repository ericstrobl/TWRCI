MCC_annotations <- function(SNPs_est,SNPs_GT){
  
  MCCs = c()
  L = length(SNPs_GT)
  
  # all_SNPs = setdiff(unlist(SNPs_GT), SNPs_GT[[L]])
  all_SNPs = unlist(SNPs_GT)
  
  if (length(SNPs_est) < length(SNPs_GT)){
    SNPs_est[length(SNPs_GT)] = list(NULL)
  }

  for (s in 1:length(SNPs_GT)){
    # SNP = setdiff(SNPs_GT[[s]], SNPs_GT[[L]])
    # an = setdiff(SNPs_est[[s]], SNPs_GT[[L]])
    SNP = SNPs_GT[[s]]
    an = SNPs_est[[s]]
    
    TP = length(intersect(an,SNP))
    FP = length(setdiff(an,SNP))
    
    TN = length(intersect(setdiff(all_SNPs,an),setdiff(all_SNPs,SNP)))
    FN = length(intersect(setdiff(all_SNPs,an),SNP))

    MCCd = sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
    if (MCCd == 0){
      MCCd = 1
    } 
    
    MCC = (TP*TN - FP*FN)/MCCd
    
    MCCs = c(MCCs, MCC)
  }
  
  return(mean(MCCs))
  
}

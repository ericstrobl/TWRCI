graph_MCC <- function(G,G_est){
  
  p = ncol(G)
  
  
  edge = which(G)
  edge_est = which(G_est)
  
  TP = length(intersect(edge_est,edge))
  FP = length(setdiff(edge_est,edge))
  
  TN = length(intersect(setdiff(1:p^2,edge_est),setdiff(1:p^2,edge)))
  FN = length(intersect(setdiff(1:p^2,edge_est),edge))
    
  MCCd = sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if (MCCd == 0){
    MCCd = 1
  } 
  
  MCC = (TP*TN - FP*FN)/MCCd
    
  return(MCC)
  
}
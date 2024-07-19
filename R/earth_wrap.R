earth_wrap <- function(x,y,z,SNPs=NULL,suffStat){
  
  x1=suffStat$data[,x];
  y1=suffStat$data[,y];
  if (is.null(SNPs)){
    z1=cbind(suffStat$data[,z,drop=FALSE],suffStat$batches)
  } else{
    z1=cbind(suffStat$data[,z],suffStat$SNP_data[,SNPs],suffStat$batches)
  }
  
  out = earth_test(x1,y1,z1)
  
  return(out$p)
  
}
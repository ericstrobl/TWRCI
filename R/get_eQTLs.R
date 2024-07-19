get_eQTLs <- function(gene_data, SNP_data, ng){
  
  p = ncol(gene_data)
  eQTLs = vector("list",p)
  
  cors = abs(cor(SNP_data,gene_data))
  
  for (g in 1:p){
    eQTLs[[g]] = which(cors[,g] > apply(cors[,-g], 1, max))
    snp = intersect(eQTLs[[g]],ng[[g]])
    if (length(snp)==0){
      eQTLs[g] = list(NULL)
    } else{
      eQTLs[[g]] = snp
    }
  }
  
  return(eQTLs)
  
}
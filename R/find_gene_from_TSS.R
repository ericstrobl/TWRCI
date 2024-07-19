find_gene_from_TSS <- function(query1,gene_names,gene_info,dis=Inf){
  # TSS contains 3 columns: chromosome, start sites, ENSG IDs
  # query1 contains 2 columns: chromosome, position
  
  require(FNN)
  require(biomaRt)
  
  # TSS = read.csv("hg19_ensembl_TSS.csv")
  # ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  # gene_info <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), 
  #                    filters ='ensembl_gene_id', values = gene_names, mart = ensembl)
  
  TSS = cbind(gene_info$chromosome_name,gene_info$start_position,gene_info$ensembl_gene_id)
  
  iT = which(TSS[,3] %in% gene_names)
  TSS = TSS[iT,]
  
  nearest_gene = rep(0,nrow(query1))
  for (c in 1:length(unique(query1[,1]))){
    ic = which(query1[,1] == c)
    icT = which(TSS[,1]==c)
    if (length(icT)>0){
      ig = knnx.index(as.numeric(TSS[icT,2]),as.numeric(query1[ic,2]),k=1)
      ig1 = which(abs(as.numeric(TSS[icT[ig],2])-as.numeric(query1[ic,2]))<dis)
      nearest_gene[ic[ig1]] = TSS[icT[ig[ig1]],3]
    }
    # print(abs(as.numeric(TSS[icT[ig],2])-as.numeric(query1[ic,2])))
    # ig = ig[which(abs(as.numeric(TSS[icT[ig],2])-as.numeric(query1[ic,2]))<dis)]
  }

  SNPs = vector("list",length(gene_names)+1)
  for (g in 1:length(gene_names)){
    SNPs[[g]] = which(nearest_gene==gene_names[g])
  }
  
  return(list(genes = gene_names,SNPs = SNPs))
  
}
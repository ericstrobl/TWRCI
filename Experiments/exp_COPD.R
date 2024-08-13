require(pcalg)
require(ieugwasr)
require(qvalue)

reps = 10

Gs = vector("list",reps)
TWRCI_res = Gs
locs_res = Gs
cis_res = Gs
ctwas_res = Gs
cis_eQTLs_res = Gs
coloc_ABF_res = Gs
coloc_susie_res = Gs

SIGNET_res = Gs
RCI_res = Gs
GRCI_res = Gs
CausalCell_res = Gs

load("COPD_exp_SNP_data_controlled2.RData")
gene_names =gsub("\\..*","",colnames(ge_f))
colnames(ge_f) = gene_names

# keep = remove_multicollinear(SNP_data_f)
out = get_genes_split(SNP_data_f,cbind(ge_f,target_f),leads_f$p)
SNP_data_f = SNP_data_f[,out$ix,drop=FALSE]
leads_f = leads_f[out$ix,,drop=FALSE]

p = ncol(ge_f)+1
ge_fs = ge_f[,setdiff(out$i1,p)]
p = ncol(ge_fs)+1

nuisance_fs = ge_f[,setdiff(out$i2,p)]

# stratified cross-validation indices
ncv = 10
folds = rep(0,length.out=length(target_f))
folds[target_f>=0.5] = rep(1:ncv,length.out=sum(target_f>=0.5))
folds[target_f<0.5] = rep(1:ncv,length.out=sum(target_f<0.5))

#TSS = read.csv("hg19_ensembl_TSS.csv")
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
#gene_info <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), 
#                    filters ='ensembl_gene_id', values = colnames(ge_fs), mart = ensembl)
#save(file="gene_info_COPD.RData",gene_info)
load("gene_info_COPD.RData")

# get descendants for silver standard
bb = kegg.gsets(species = "hsa", id.type="entrez", check.new=FALSE) 
library("AnnotationDbi")
library("org.Hs.eg.db")
EnsID = mapIds(org.Hs.eg.db,
               keys=bb$kg.sets$`hsa04660 T cell receptor signaling pathway`,
               column="ENSEMBL",
               keytype="ENTREZID",
               multiVals="first")
desc = which(colnames(ge_f) %in% EnsID)
ps= c()
for (d in desc){
  ps = c(ps, cor.test(ge_f[,d],target_f)$p.value)
}
desc = desc[which(ps<0.05)]

# run stratified cross-validation
for (i in 1:ncv){
  print(i)
  
  ge_fst = ge_fs[folds!=i,]
  SNP_data_ft = SNP_data_f[folds!=i,]
  target_ft = target_f[folds!=i]
  age_ft = age_f[folds!=i]
  nuisance_ft = nuisance_fs[folds!=i,]
  
  # betas
  prod = SNP_data_ft * target_ft
  SY = colSums(prod)
  corZY = colMeans(prod)
  corZY_var = apply(prod,2,var)/nrow(SNP_data_ft)
  batches = age_ft
  
  ### annotation
  ptm <- proc.time()
  aa = TWRCI(ge_fst,SNP_data_ft,target_ft,nuisance_ft,batches)
  TWRCI_res[[i]]$time = (proc.time() - ptm)[3]
  TWRCI_res[[i]]$annot = aa
  # aa = TWRCI_res[[i]]$annot
  
  ### locations
  query1 = cbind(leads_f$chr,leads_f$position)
  
  ### TSS
  ptm <- proc.time()
  ng_1 = find_gene_from_TSS(query1,colnames(ge_fs),gene_info)
  locs_res[[i]]$time = (proc.time() - ptm)[3]
  locs_res[[i]]$annot = ng_1

  ## cis-SNPs
  ptm <- proc.time()
  ng = find_gene_from_TSS(query1,colnames(ge_fs),gene_info,500000)
  cis_res[[i]]$time = (proc.time() - ptm)[3]
  cis_res[[i]]$annot = ng

  ### cTWAS
  ptm <- proc.time()
  cTWAS = ctwas(SNP_data_ft,ge_fst,target_ft,ng$SNPs,cbind(nuisance_ft,age_ft)) # directories in this function need to be customized by following https://xinhe-lab.github.io/ctwas/articles/ctwas_summary_statistics.html
  ctwas_res[[i]]$time = (proc.time() - ptm)[3]
  ctwas_res[[i]]$annot = cTWAS
  
  ### cis-eQTLs
  ptm <- proc.time()
  cis_eQTLs = get_eQTLs(cbind(ge_fst,target_ft),SNP_data_ft,ng$SNPs) # correlated and nearby TSS
  cis_eQTLs_res[[i]]$time = (proc.time() - ptm)[3]
  cis_eQTLs_res[[i]]$annot = cis_eQTLs
  
  ### colocalization, ABF
  ptm <- proc.time()
  CL = get_coloc_ss(ge_fst,SNP_data_ft,batches,corZY,corZY_var,ng_1$SNPs)
  coloc_ABF_res[[i]]$time = (proc.time() - ptm)[3]
  coloc_ABF_res[[i]]$annot = CL
  
  ### colocalization, SuSIE
  ptm <- proc.time()
  CL_susie = get_coloc_ss(ge_fst,SNP_data_ft,batches,corZY,corZY_var,ng_1$SNPs,method="susie")
  coloc_susie_res[[i]]$time = (proc.time() - ptm)[3]
  coloc_susie_res[[i]]$annot = CL_susie

  
  ### network reconstruction
  suffStat = list(); suffStat$batches = normalizeData(cbind(nuisance_ft,age_ft));
  suffStat$data = normalizeData(cbind(ge_fst,target_ft));
  suffStat$SNP_data = normalizeData(SNP_data_ft)
  
  # TWRCI
  print("TWRCI_graph")
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, SNPs=aa$SNPs, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  TWRCI_res[[i]]$skel = skel
  TWRCI_res[[i]]$time_skel = (proc.time() - ptm)[3]
  # skel = TWRCI_res[[i]]$skel
  
  ptm <- proc.time()
  G_est =  TWRCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[aa$K[k],aa$K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=aa$SNPs, alpha=0.05, root=FALSE) # check graph consistency with binary target
  # colnames(G_est1) = gsub("\\..*","",colnames(ge_fs))
  TWRCI_res[[i]]$time_graph = (proc.time() - ptm)[3]
  TWRCI_res[[i]]$G_est = G_est
  # G_est = TWRCI_res[[i]]$G_est
  
  Xte = SNP_data_f[folds==i,]
  Yte = cbind(ge_fs[folds==i,],target_f[folds==i,])
  Bte = cbind(nuisance_fs[folds==i,],age_f[folds==i])
  batchest = cbind(nuisance_ft,age_ft)
  
  TWRCI_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,aa$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,ng_1$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,ng$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,cTWAS,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,cis_eQTLs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,CL,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,CL_susie,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)

  TWRCI_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,aa$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  locs_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,ng_1$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,ng$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  ctwas_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,cTWAS,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,cis_eQTLs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,CL,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_susie_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,CL_susie,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])

  TWRCI_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],aa$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  locs_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng_1$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  ctwas_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cTWAS,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_silver =  MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cis_eQTLs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_susie_res[[i]]$MACR_silver = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL_susie,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  
  
  ## get graph for SIGNET below
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, SNPs = ng_1$SNPs, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  SIGNET_res[[i]]$skel = skel
  SIGNET_res[[i]]$time_skel = (proc.time() - ptm)[3]
  # skel = SIGNET_res[[i]]$skel
  
  # SIGNET
  print("SIGNET_graph")
  ptm <- proc.time()
  G_est = SIGNET(ge_fst,batches,SNP_data_ft,ng_1$SNPs)>0
  K = c(extract_causal_order(G_est),p)
  G_est =  SIGNET_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=ng_1$SNPs, alpha=0.05, root=FALSE) # check graph consistency with binary target
  SIGNET_res[[i]]$time_graph = (proc.time() - ptm)[3]
  SIGNET_res[[i]]$G_est = G_est
  # G_est = SIGNET_res[[i]]$G_est
  
  TWRCI_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,aa$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,ng_1$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,ng$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,cTWAS,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,cis_eQTLs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,CL,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,CL_susie,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,aa$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  locs_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,ng_1$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,ng$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  ctwas_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,cTWAS,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,cis_eQTLs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,CL,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_susie_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,CL_susie,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])

  TWRCI_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],aa$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  locs_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng_1$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  ctwas_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cTWAS,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_silver_SIGNET =  MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cis_eQTLs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_susie_res[[i]]$MACR_silver_SIGNET = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL_susie,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  
  
  ## get graph for RCI, GRCI and PC below
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  RCI_res[[i]]$skel = skel
  RCI_res[[i]]$time_skel = (proc.time() - ptm)[3]
  # skel = RCI_res[[i]]$skel

  ### RCI
  print("RCI_graph")
  ptm <- proc.time()
  K = c(RCI_orig(ge_fst),p)
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05)
  # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  RCI_res[[i]]$time_order = (proc.time() - ptm)[3]
  RCI_res[[i]]$G_est = G_est
  # G_est = RCI_res[[i]]$G_est
  
  TWRCI_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,aa$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,ng_1$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,ng$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,cTWAS,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,cis_eQTLs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,CL,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,CL_susie,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,aa$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  locs_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,ng_1$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,ng$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  ctwas_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,cTWAS,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,cis_eQTLs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,CL,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_susie_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,CL_susie,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])

  TWRCI_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],aa$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  locs_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng_1$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  ctwas_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cTWAS,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_silver_RCI =  MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cis_eQTLs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_susie_res[[i]]$MACR_silver_RCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL_susie,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  
  ### GRCI
  print("GRCI_graph")
  ptm <- proc.time()
  K = c(GRCI_ANM(ge_fst),p)
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05,root=FALSE)
  # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  GRCI_res[[i]]$time_order = (proc.time() - ptm)[3]
  GRCI_res[[i]]$G_est = G_est
  # G_est = GRCI_res[[i]]$G_est

  TWRCI_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,aa$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,ng_1$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,ng$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,cTWAS,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,cis_eQTLs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,CL,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,CL_susie,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)

  TWRCI_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,aa$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  locs_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,ng_1$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,ng$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  ctwas_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,cTWAS,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,cis_eQTLs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,CL,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_susie_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,CL_susie,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])

  TWRCI_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],aa$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  locs_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng_1$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  ctwas_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cTWAS,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_silver_GRCI =  MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cis_eQTLs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_susie_res[[i]]$MACR_silver_GRCI = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL_susie,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  

  ### CausalCell / PC
  print("PC_graph")
  ptm <- proc.time()
  G_orient = my_udag2pdag(RCI_res[[i]]$skel$G,RCI_res[[i]]$skel$sepset)
  K = c(extract_causal_order(G_orient[-p,-p]),p) # rerun with causal order to control for level of sparsity before computing SHD
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05,root=FALSE)
  # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  CausalCell_res[[i]]$time_order = (proc.time() - ptm)[3]
  CausalCell_res[[i]]$G_est = G_est
  # G_est = CausalCell_res[[i]]$G_est

  TWRCI_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,aa$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,ng_1$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,ng$SNPs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,cTWAS,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,cis_eQTLs,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,CL,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,CL_susie,SNP_data_ft,cbind(ge_fst,target_ft),batchest,Xte,Yte,Bte)

  TWRCI_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,aa$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  locs_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,ng_1$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,ng$SNPs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  ctwas_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,cTWAS,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,cis_eQTLs,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,CL,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])
  coloc_susie_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,CL_susie,SNP_data_ft,ge_fst,batchest,target_ft,SNP_data_f[folds==i,],ge_f[folds==i,],Bte,target_f[folds==i,])

  TWRCI_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],aa$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  locs_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng_1$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],ng$SNPs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  ctwas_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cTWAS,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  cis_eQTLs_res[[i]]$MACR_silver_CC =  MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],cis_eQTLs,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_ABF_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  coloc_susie_res[[i]]$MACR_silver_CC = MACR_silver(target_f[folds==i],G_est,ge_fs[folds==i,],SNP_data_f[folds==i,],CL_susie,cbind(nuisance_fs,age_f)[folds==i,],desc,ge_f[folds==i,])
  
  
  save(file="Results_COPD.RData", Gs, TWRCI_res, locs_res, cis_res, ctwas_res, cis_eQTLs_res, coloc_ABF_res, coloc_susie_res, 
       SIGNET_res, RCI_res, GRCI_res, CausalCell_res)
}


# annotation
imax = 10
res_mat = array(0,c(imax,7,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MACR_annot_graph
  res_mat[i,2,1] = locs_res[[i]]$MACR_annot_graph
  res_mat[i,3,1] = cis_res[[i]]$MACR_annot_graph
  res_mat[i,4,1] = ctwas_res[[i]]$MACR_annot_graph
  res_mat[i,5,1] = cis_eQTLs_res[[i]]$MACR_annot_graph
  res_mat[i,6,1] = coloc_ABF_res[[i]]$MACR_annot_graph
  res_mat[i,7,1] = coloc_susie_res[[i]]$MACR_annot_graph
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,4,2] = ctwas_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,5,2] = cis_eQTLs_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,6,2] = coloc_ABF_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,7,2] = coloc_susie_res[[i]]$MACR_annot_graph_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,2,3] = locs_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,3,3] = cis_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,4,3] = ctwas_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,5,3] = cis_eQTLs_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,6,3] = coloc_ABF_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,7,3] = coloc_susie_res[[i]]$MACR_annot_graph_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,4,4] = ctwas_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,5,4] = cis_eQTLs_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,6,4] = coloc_ABF_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,7,4] = coloc_susie_res[[i]]$MACR_annot_graph_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MACR_annot_graph_CC
  res_mat[i,2,5] = locs_res[[i]]$MACR_annot_graph_CC
  res_mat[i,3,5] = cis_res[[i]]$MACR_annot_graph_CC
  res_mat[i,4,5] = ctwas_res[[i]]$MACR_annot_graph_CC
  res_mat[i,5,5] = cis_eQTLs_res[[i]]$MACR_annot_graph_CC
  res_mat[i,6,5] = coloc_ABF_res[[i]]$MACR_annot_graph_CC
  res_mat[i,7,5] = coloc_susie_res[[i]]$MACR_annot_graph_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1*t(apply(res_mat,c(2,3),sd))/sqrt(10)
t(apply(res_mat,c(2,3),mean))-1*t(apply(res_mat,c(2,3),sd))/sqrt(10)


# CRCE
imax =10
res_mat = array(0,c(imax,7,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MACR_CRCE
  res_mat[i,2,1] = locs_res[[i]]$MACR_CRCE
  res_mat[i,3,1] = cis_res[[i]]$MACR_CRCE
  res_mat[i,4,1] = ctwas_res[[i]]$MACR_CRCE
  res_mat[i,5,1] = cis_eQTLs_res[[i]]$MACR_CRCE
  res_mat[i,6,1] = coloc_ABF_res[[i]]$MACR_CRCE
  res_mat[i,7,1] = coloc_susie_res[[i]]$MACR_CRCE
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,4,2] = ctwas_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,5,2] = cis_eQTLs_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,6,2] = coloc_ABF_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,7,2] = coloc_susie_res[[i]]$MACR_CRCE_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MACR_CRCE_RCI
  res_mat[i,2,3] = locs_res[[i]]$MACR_CRCE_RCI
  res_mat[i,3,3] = cis_res[[i]]$MACR_CRCE_RCI
  res_mat[i,4,3] = ctwas_res[[i]]$MACR_CRCE_RCI
  res_mat[i,5,3] = cis_eQTLs_res[[i]]$MACR_CRCE_RCI
  res_mat[i,6,3] = coloc_ABF_res[[i]]$MACR_CRCE_RCI
  res_mat[i,7,3] = coloc_susie_res[[i]]$MACR_CRCE_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,4,4] = ctwas_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,5,4] = cis_eQTLs_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,6,4] = coloc_ABF_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,7,4] = coloc_susie_res[[i]]$MACR_CRCE_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MACR_CRCE_CC
  res_mat[i,2,5] = locs_res[[i]]$MACR_CRCE_CC
  res_mat[i,3,5] = cis_res[[i]]$MACR_CRCE_CC
  res_mat[i,4,5] = ctwas_res[[i]]$MACR_CRCE_CC
  res_mat[i,5,5] = cis_eQTLs_res[[i]]$MACR_CRCE_CC
  res_mat[i,6,5] = coloc_ABF_res[[i]]$MACR_CRCE_CC
  res_mat[i,7,5] = coloc_susie_res[[i]]$MACR_CRCE_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1*t(apply(res_mat,c(2,3),sd))/sqrt(10)
t(apply(res_mat,c(2,3),mean))-1*t(apply(res_mat,c(2,3),sd))/sqrt(10)


# silver standard
imax =10
res_mat = array(0,c(imax,7,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MACR_silver
  res_mat[i,2,1] = locs_res[[i]]$MACR_silver
  res_mat[i,3,1] = cis_res[[i]]$MACR_silver
  res_mat[i,4,1] = ctwas_res[[i]]$MACR_silver
  res_mat[i,5,1] = cis_eQTLs_res[[i]]$MACR_silver
  res_mat[i,6,1] = coloc_ABF_res[[i]]$MACR_silver
  res_mat[i,7,1] = coloc_susie_res[[i]]$MACR_silver
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MACR_silver_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MACR_silver_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MACR_silver_SIGNET
  res_mat[i,4,2] = ctwas_res[[i]]$MACR_silver_SIGNET
  res_mat[i,5,2] = cis_eQTLs_res[[i]]$MACR_silver_SIGNET
  res_mat[i,6,2] = coloc_ABF_res[[i]]$MACR_silver_SIGNET
  res_mat[i,7,2] = coloc_susie_res[[i]]$MACR_silver_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MACR_silver_RCI
  res_mat[i,2,3] = locs_res[[i]]$MACR_silver_RCI
  res_mat[i,3,3] = cis_res[[i]]$MACR_silver_RCI
  res_mat[i,4,3] = ctwas_res[[i]]$MACR_silver_RCI
  res_mat[i,5,3] = cis_eQTLs_res[[i]]$MACR_silver_RCI
  res_mat[i,6,3] = coloc_ABF_res[[i]]$MACR_silver_RCI
  res_mat[i,7,3] = coloc_susie_res[[i]]$MACR_silver_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MACR_silver_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MACR_silver_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MACR_silver_GRCI
  res_mat[i,4,4] = ctwas_res[[i]]$MACR_silver_GRCI
  res_mat[i,5,4] = cis_eQTLs_res[[i]]$MACR_silver_GRCI
  res_mat[i,6,4] = coloc_ABF_res[[i]]$MACR_silver_GRCI
  res_mat[i,7,4] = coloc_susie_res[[i]]$MACR_silver_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MACR_silver_CC
  res_mat[i,2,5] = locs_res[[i]]$MACR_silver_CC
  res_mat[i,3,5] = cis_res[[i]]$MACR_silver_CC
  res_mat[i,4,5] = ctwas_res[[i]]$MACR_silver_CC
  res_mat[i,5,5] = cis_eQTLs_res[[i]]$MACR_silver_CC
  res_mat[i,6,5] = coloc_ABF_res[[i]]$MACR_silver_CC
  res_mat[i,7,5] = coloc_susie_res[[i]]$MACR_silver_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1*t(apply(res_mat,c(2,3),sd))/sqrt(10)
t(apply(res_mat,c(2,3),mean))-1*t(apply(res_mat,c(2,3),sd))/sqrt(10)

# timing annotation
imax = 10
res_mat = matrix(0,imax,7)
for (i in 1:imax){
  print(i)
  res_mat[i,1] = TWRCI_res[[i]]$time
  res_mat[i,2] = locs_res[[i]]$time
  res_mat[i,3] = cis_res[[i]]$time
  res_mat[i,4] = ctwas_res[[i]]$time
  res_mat[i,5] = cis_eQTLs_res[[i]]$time
  res_mat[i,6] = coloc_ABF_res[[i]]$time
  res_mat[i,7] = coloc_susie_res[[i]]$time
}
print(colMeans(res_mat))

# timing graph
imax = 10
res_mat = matrix(0,imax,5)
for (i in 1:imax){
  print(i)
  res_mat[i,1] = TWRCI_res[[i]]$time_skel + TWRCI_res[[i]]$time_graph
  res_mat[i,2] = SIGNET_res[[i]]$time_skel + SIGNET_res[[i]]$time_graph
  res_mat[i,3] = RCI_res[[i]]$time_skel + RCI_res[[i]]$time_order
  res_mat[i,4] = RCI_res[[i]]$time_skel + GRCI_res[[i]]$time_order
  res_mat[i,5] = RCI_res[[i]]$time_skel + CausalCell_res[[i]]$time_order
}
print(colMeans(res_mat))


### run all samples
aa = TWRCI(ge_fs,SNP_data_f,target_f,nuisance_fs,age_f)

# causal graph
suffStat = list(); suffStat$batches = normalizeData(cbind(nuisance_fs,age_f)); 
suffStat$data = normalizeData(cbind(ge_fs,target_f)); 
suffStat$SNP_data = normalizeData(SNP_data_f)

skel = my_skeleton_p(suffStat, p, SNPs=aa$SNPs, alpha=0.05)
G_est = skel$G
for (k in seq_len(p)){
  G_est[aa$K[k],aa$K[seq_len(k-1)]]=FALSE
}
G_est = screen_anc_Y(suffStat, G_est, SNPs=aa$SNPs,alpha=0.05,root=FALSE)
plot(as(G_est,"graphNEL"))

# pie chart, proportions
iS = isAncAll(G_est,1:p,p)
props = c()
for (k in iS){
  props = c(props, length(aa$SNPs[[k]])/length(unlist(aa$SNPs[iS])))
}
print(props)

# histogram of distances
iS = isAncAll(G_est,1:(p-1),p)
gene_names = colnames(ge_fs)[iS]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
gene_info <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), 
                   filters ='ensembl_gene_id', values = gene_names, mart = ensembl)

TSS = cbind(gene_info$chromosome_name,gene_info$start_position,gene_info$ensembl_gene_id)

iT = which(TSS[,3] %in% gene_names)
TSS = TSS[iT,]

alphas = rep(0,ncol(SNP_data_f))
mmod = CV_LRR2_cor(SNP_data_f[,unlist(aa$SNPs[iS])],target_f)
alphas[unlist(aa$SNPs[iS])] = mmod$betas

dis = rep(0,ncol(SNP_data_f))
genes = setdiff(iS,p)
for (s in 1:length(genes)){
  chr = leads_f$chr[aa$SNPs[[genes[s]]]]
  icT = which(chr!= TSS[s,1])
  dis[aa$SNPs[[genes[s]]][icT]] = Inf
  icT = which(chr== TSS[s,1])
  dis[aa$SNPs[[genes[s]]][icT]] = abs(as.numeric(TSS[s,2])-as.numeric(leads_f$position[aa$SNPs[[genes[s]]][icT]]))
}
dis = dis[unlist(aa$SNPs[genes])]


dist = classInt::classIntervals(dis[dis!=Inf], 50,style="quantile")$brks
perc = c()
for (d in dist){
  perc = c(perc,(1-(sum(dis<d)/length(dis)))*100)
}
print(cbind(dist,perc))


# moving average
iS = union(p,iS)
alphas = CV_LRR2(SNP_data_f[,unlist(aa$SNPs[iS])],target_f)$betas
alphas = alphas[-(1:length(aa$SNPs[p]))] # remove alphas for Y
A = dis[dis!=Inf][order(dis[dis!=Inf])]
B = ma(alphas[dis!=Inf][order(dis[dis!=Inf])],n=10)
print(cbind(A,B))

print(mean(alphas[dis==Inf]))


# clustering
require(uwot)

#load("COPD_UMAP_age_all.RData")

iS = isAncAll(G_est,1:(p-1),p)
ge_fs_age = normalizeData(earth(age_f,ge_fs)$residuals)

data.umap = normalizeData(uwot::umap(ge_fs_age[,setdiff(iS,p)],y=factor(target_f>1)))
plot(data.umap[,1],data.umap[,2],col=(target_f>1)+1)
require(fastcluster)
cl <- hclust.vector(data.umap[target_f>1,], method="ward")
plot(rev(cl$height)[1:10])
k = pathviewr:::find_curve_elbow(cbind(1:10,rev(cl$height)[1:10]))
ix = cutree(cl, k = k)
plot(data.umap[target_f>1,],col=ix)
colMeans(mm$CRCE[target_f>1,][ix==1,])
data.umap[target_f>1,][ix==1,]
data.umap[target_f>1,][ix==2,]
data.umap[target_f>1,][ix==3,]
data.umap[target_f>1,][ix==4,]

# scaling
iS = isAncAll(G_est,1:(p-1),p)
suffStat$batches = normalizeData(cbind(nuisance_fs,age_f))
mm = estimate_CRCEs(ge_fs,age_f,G_est,aa$SNPs,SNP_data_f,target_f,iS)
colMeans(mm$CRCE[target_f>1,])/mean(abs(colMeans(mm$CRCE[target_f>1,])))
colMeans(mm$CRCE[target_f>1,][ix==1,])/mean(abs(colMeans(mm$CRCE[target_f>1, ][ix==1,])))
colMeans(mm$CRCE[target_f>1,][ix==2,])/mean(abs(colMeans(mm$CRCE[target_f>1, ][ix==2,])))
colMeans(mm$CRCE[target_f>1,][ix==3,])/mean(abs(colMeans(mm$CRCE[target_f>1, ][ix==3,])))
colMeans(mm$CRCE[target_f>1,iS][ix==4,])/mean(abs(colMeans(mm$CRCE[target_f>1, ][ix==4,])))

# gene expression levels of clusters
for(k in 1:4){
  mean = colMeans(ge_fs[target_f>1,iS][ix==k,])-colMeans(ge_fs[target_f<=1,iS])
  sd = sqrt(apply(ge_fs[target_f>1,iS][ix==k,],2,var)/sum(ix==k) + apply(ge_fs[target_f<=1,iS],2,var)/sum(target_f<=1))
  print(mean)
  print(mean+1.96*sd)
  print(mean-1.96*sd)
  print("--------")
}

# histogram of correlations, confounding
rs = cor(SNP_data_f[,aa$SNPs[[p]]],ge_fs)
ts=rs*sqrt((nrow(ge_fs)-2)/(1-rs^2))
hist( (1-pt(abs(ts),nrow(ge_fs)-2))*2)
hh = hist( (1-pt(abs(ts),nrow(ge_fs)-2))*2)
hh$mids
hh$counts
ks.test((1-pt(abs(ts),nrow(ge_fs)-2))*2,"punif",0,1)

# trans-eQTLs
require(FNN)
trans = read.csv("trans_eQTLs.csv") # https://eqtlgen.org/trans-eqtls.html
anc = isAncAll(G_est,1:(p-1),p)
leads_fa = leads_f[unlist(aa$SNPs[anc]),]

dists = c()
for (c in 1:22){
  it = which(trans$SNPChr == c)
  il = which(leads_fa$chr == as.character(c))
  
  pos = c(trans$SNPPos[it],leads_fa$position[il])
  trans$SNPPos[it] = trans$SNPPos[it]/sd(pos)
  leads_fa$position[il] = leads_fa$position[il]/sd(pos)
}

dists = c()
for (c in 1:22){
  it = which(trans$SNPChr == c)
  il = which(leads_fa$chr == as.character(c))
  
  dist = knnx.dist(trans$SNPPos[it],leads_fa$position[il],k=1)
  dists = c(dists,dist)
}
dist_stat = median(dists)

dist_median = c()
for (p in 1:10000){
  print(p)
  r = length(trans$SNPPos)
  transt = trans
  transt$SNPPos = transt$SNPPos[sample(1:r,r,replace=FALSE)]
  dists = c()
  for (c in 1:22){
    it = which(transt$SNPChr == c)
    il = which(leads_fa$chr == as.character(c))
    
    dist = knnx.dist(transt$SNPPos[it],leads_fa$position[il],k=1)
    dists = c(dists,dist)
  }
  dist_median = c(dist_median, median(dists))
}
mean(dist_median<dist_stat)
mean(dist_median/dist_stat)
mean(dist_median/dist_stat) + 1.96*sd(dist_median/dist_stat)/sqrt(10000)
mean(dist_median/dist_stat) - 1.96*sd(dist_median/dist_stat)/sqrt(10000)

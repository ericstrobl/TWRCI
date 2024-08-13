require(pcalg)
require(ieugwasr)
require(qvalue)

dzs = c("ieu-b-5067","ieu-b-4967","ieu-b-4972","ieu-b-4971","ieu-b-4973","ieu-a-1187",
        "ieu-b-4965","ieu-b-5063","ieu-b-4956","ieu-b-4953")

reps = 500

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

nsamps = 500
var_file = rep(1:10,length.out=reps)
for (i in 51:100){
  print(i)
  load(paste(dzs[var_file[i]],"_preprocessed.RData",sep=""))
  
  ib = sample(1:nrow(data),nsamps,replace=TRUE)
  keep1 = which(apply(data[ib,],2,sd)>1E-20)
  print(length(keep1))
  keep = remove_multicollinear(data[ib,keep1])

  datat = data[ib,keep1[keep]]
  

  p = 30
  SNPs = generate_SNP_list(p,ncol(datat),leads) # assign SNPs to each gene
  G = generate_DAG(p,SNPs$SNPs,bulk_nbatch=1) # generate DAG over RNA-seq, ancestral to Y
  Gs[[i]] = G
  plot(as(as.matrix(G$graph),"graphNEL"))
  RNAseq = sample_data(G,SNPs$SNPs,datat,p=p,nsamps=nsamps) # generate RNA-seq data
  
  # betas
  prod = normalizeData(RNAseq$SNP_data) * c(normalizeData(RNAseq$Y))
  SY = colSums(prod)
  corZY = colMeans(prod)
  corZY_var = apply(prod,2,var)/nsamps
  LD = cor(normalizeData(RNAseq$SNP_data))
  
  
  ### annotation
  ptm <- proc.time()
  aa = TWRCI(RNAseq$X,RNAseq$SNP_data,RNAseq$Y,RNAseq$batch)
  TWRCI_res[[i]]$time = (proc.time() - ptm)[3]
  # aa = RCI11_Y(RNAseq$X,RNAseq$Y,RNAseq$SNP_data,RNAseq$batch)
  TWRCI_res[[i]]$MCC = MCC_annotations(aa$SNPs,SNPs$SNPs)
  TWRCI_res[[i]]$rank = rank_annot(aa$SNPs,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G,K=aa$K)
  TWRCI_res[[i]]$MCC_HP = MCC_annotations(aa$SNPs,SNPs$SNPs,p)
  TWRCI_res[[i]]$rank_HP = rank_annot(aa$SNPs,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G,p)
  
  ### locations
  # time not applicable here, time = 0 automatically
  ng_1 = get_cis_SNPs_window(SNPs$locs[-p],leads[keep1[keep],],window=Inf)
  ng_1[p] = list(NULL)
  locs_res[[i]]$MCC = MCC_annotations(ng_1,SNPs$SNPs)
  locs_res[[i]]$rank = rank_annot(ng_1,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)
  
  ### cis-SNPs
  ptm <- proc.time()
  cis = get_cis_SNPs_window(SNPs$locs[-p],leads[keep1[keep],]) # cis SNPs among SNPs passing 5e-5
  cis[p] = list(NULL)
  cis_res[[i]]$time = (proc.time() - ptm)[3]
  cis_res[[i]]$MCC = MCC_annotations(cis,SNPs$SNPs)
  cis_res[[i]]$rank = rank_annot(cis,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)

  ### cTWAS
  ptm <- proc.time()
  cTWAS = ctwas(RNAseq$SNP_data,RNAseq$X,RNAseq$Y,cis,RNAseq$batch)
  ctwas_res[[i]]$time = (proc.time() - ptm)[3]
  ctwas_res[[i]]$MCC = MCC_annotations(cTWAS,SNPs$SNPs)
  ctwas_res[[i]]$rank = rank_annot(cTWAS,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)
  ctwas_res[[i]]$MCC_HP = MCC_annotations(cTWAS,SNPs$SNPs,p)
  ctwas_res[[i]]$rank_HP = rank_annot(cTWAS,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G,p)
  
  ### eQTLs
  ptm <- proc.time()
  cis_eQTLs = get_eQTLs(cbind(RNAseq$X,RNAseq$Y),RNAseq$SNP_data,cis)
  cis_eQTLs_res[[i]]$time = (proc.time() - ptm)[3]
  cis_eQTLs_res[[i]]$MCC = MCC_annotations(cis_eQTLs,SNPs$SNPs)
  cis_eQTLs_res[[i]]$rank = rank_annot(cis_eQTLs,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)
  
  ## colocalization, ABF
  ptm <- proc.time()
  CL = get_coloc_ss(RNAseq$X,RNAseq$SNP_data,RNAseq$batch,corZY,corZY_var,ng_1) # Predixcan, MV-IWAS, 2SLS
  coloc_ABF_res[[i]]$time = (proc.time() - ptm)[3]  
  coloc_ABF_res[[i]]$MCC = MCC_annotations(CL,SNPs$SNPs)
  coloc_ABF_res[[i]]$rank = rank_annot(CL,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)
  
  ## colocalization, SuSIE
  ptm <- proc.time()
  CL_susie = get_coloc_ss(RNAseq$X,RNAseq$SNP_data,RNAseq$batch,corZY,corZY_var,ng_1,method="susie") # Predixcan, MV-IWAS, 2SLS
  coloc_susie_res[[i]]$time = (proc.time() - ptm)[3]
  coloc_susie_res[[i]]$MCC = MCC_annotations(CL_susie,SNPs$SNPs) ##
  coloc_susie_res[[i]]$rank = rank_annot(CL_susie,SNPs$SNPs,RNAseq$SNP_data,cbind(RNAseq$X,RNAseq$Y),G)
  
  
  ### network reconstruction
  suffStat = list(); suffStat$data = normalizeData(cbind(RNAseq$X,RNAseq$Y)); suffStat$SNP_data = normalizeData(RNAseq$SNP_data)
  suffStat$batches = rep(1,nrow(RNAseq$X))
  
  G_tr = G$graph
  G_tr[p,]=FALSE;G_tr[,p]=FALSE
  G_tr[isAncAll(G$graph,1:(p-1),p),p] = TRUE
  
  
  # TWRCI
  print("TWRCI_graph")
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, SNPs=aa$SNPs, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  TWRCI_res[[i]]$skel = skel
  TWRCI_res[[i]]$time_skel = (proc.time() - ptm)[3]
  # skel = Gs[[i]]$skel
  
  ptm <- proc.time()
  G_est =  TWRCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[aa$K[k],aa$K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=aa$SNPs, alpha=0.05) # check graph consistency with binary target
  TWRCI_res[[i]]$time_graph = (proc.time() - ptm)[3]
  TWRCI_res[[i]]$G_est = G_est
  
  TWRCI_res[[i]]$SHD = sum(G_est!=G_tr)
  TWRCI_res[[i]]$MCC_graph = graph_MCC(G_tr>0,G_est)
  
  # need test samples here
  RNAseq_te = sample_data(G,SNPs$SNPs,datat,p=p,nsamps=nsamps) # generate RNA-seq data
  
  X = normalizeData(RNAseq$SNP_data)
  Y = normalizeData(cbind(RNAseq$X,RNAseq$Y))
  B = as.matrix(RNAseq$batch)
  Xte = normalizeData(RNAseq_te$SNP_data)
  Yte = normalizeData(cbind(RNAseq_te$X,RNAseq_te$Y))
  Bte = as.matrix(RNAseq_te$batch)
  
  TWRCI_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,aa$SNPs,X,Y,B,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,ng_1,X,Y,B,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,cis,X,Y,B,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,cTWAS,X,Y,B,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,cis_eQTLs,X,Y,B,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,CL,X,Y,B,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph = MACR_annot_graph(G_est,CL_susie,X,Y,B,Xte,Yte,Bte)
  
  Yg = normalizeData(RNAseq$X)
  Yy = normalizeData(RNAseq$Y)
  Ygte = normalizeData(RNAseq_te$X)
  Yyte = normalizeData(RNAseq_te$Y)
  
  truth = MSE_CRCE(G_tr>0,SNPs$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte)
  
  TWRCI_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  locs_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_eQTLs_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_ABF_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_susie_res[[i]]$MSE_CRCE = MSE_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  
  TWRCI_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  locs_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_eQTLs_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_ABF_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_susie_res[[i]]$MACR_CRCE = MACR_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  
  # SIGNET
  ## get graph for SIGNET below
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, SNPs = ng_1, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  SIGNET_res[[i]]$skel = skel
  SIGNET_res[[i]]$time_skel = (proc.time() - ptm)[3]
  
  # SIGNET
  print("SIGNET_graph")
  ptm <- proc.time()
  G_est = SIGNET(RNAseq$X,RNAseq$batch,RNAseq$SNP_data,ng_1)>0
  K = c(extract_causal_order(G_est),p)
  G_est =  SIGNET_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=ng_1, alpha=0.05, root=FALSE) # check graph consistency with binary target
  SIGNET_res[[i]]$time_graph = (proc.time() - ptm)[3]
  SIGNET_res[[i]]$G_est = G_est
  SIGNET_res[[i]]$SHD = sum(G_est!=G$graph)
  SIGNET_res[[i]]$MCC_graph = graph_MCC(G$graph>0,G_est)
  
  TWRCI_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,aa$SNPs,X,Y,B,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,ng_1,X,Y,B,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,cis,X,Y,B,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,cTWAS,X,Y,B,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,cis_eQTLs,X,Y,B,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,CL,X,Y,B,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_SIGNET = MACR_annot_graph(G_est,CL_susie,X,Y,B,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,SNPs$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  locs_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_eQTLs_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_ABF_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_susie_res[[i]]$MSE_CRCE_SIGNET = MSE_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  
  TWRCI_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  locs_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_eQTLs_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_ABF_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_susie_res[[i]]$MACR_CRCE_SIGNET = MACR_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  
  
  ## get graph for RCI, GRCI and PC below
  ptm <- proc.time()
  skel = my_skeleton_p(suffStat, p-1, alpha=0.05)
  Gn = matrix(FALSE,p,p)
  Gn[1:(p-1),1:(p-1)] = skel$G
  skel$G = Gn
  RCI_res[[i]]$skel = skel
  RCI_res[[i]]$time_skel = (proc.time() - ptm)[3]
  
  # RCI
  print("RCI_graph")
  ptm <- proc.time()
  K = c(RCI_orig(RNAseq$X),p)
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05)
  # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  RCI_res[[i]]$time_order = (proc.time() - ptm)[3]
  RCI_res[[i]]$G_est = G_est
  
  RCI_res[[i]]$SHD = sum(G_est!=G$graph)
  RCI_res[[i]]$MCC_graph = graph_MCC(G$graph>0,G_est)
  
  TWRCI_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,aa$SNPs,X,Y,B,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,ng_1,X,Y,B,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,cis,X,Y,B,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,cTWAS,X,Y,B,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,cis_eQTLs,X,Y,B,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,CL,X,Y,B,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_RCI = MACR_annot_graph(G_est,CL_susie,X,Y,B,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  locs_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_eQTLs_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_ABF_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_susie_res[[i]]$MSE_CRCE_RCI = MSE_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  
  TWRCI_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  locs_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_eQTLs_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_ABF_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_susie_res[[i]]$MACR_CRCE_RCI = MACR_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  
  
  # GRCI
  print("GRCI_graph")
  ptm <- proc.time()
  K = c(GRCI_ANM(RNAseq$X),p)
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05)
  # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  GRCI_res[[i]]$time_order = (proc.time() - ptm)[3]
  GRCI_res[[i]]$G_est = G_est
  
  GRCI_res[[i]]$SHD = sum(G_est!=G$graph) 
  GRCI_res[[i]]$MCC_graph = graph_MCC(G$graph>0,G_est)
  
  TWRCI_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,aa$SNPs,X,Y,B,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,ng_1,X,Y,B,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,cis,X,Y,B,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,cTWAS,X,Y,B,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,cis_eQTLs,X,Y,B,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,CL,X,Y,B,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_GRCI = MACR_annot_graph(G_est,CL_susie,X,Y,B,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  locs_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_eQTLs_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_ABF_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_susie_res[[i]]$MSE_CRCE_GRCI = MSE_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  
  TWRCI_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  locs_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_eQTLs_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_ABF_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_susie_res[[i]]$MACR_CRCE_GRCI = MACR_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  
  
  # CausalCell / PC
  print("PC_graph")
  ptm <- proc.time()
  G_orient = my_udag2pdag(RCI_res[[i]]$skel$G,RCI_res[[i]]$skel$sepset)
  K = c(extract_causal_order(G_orient[-p,-p]),p) # rerun with causal order to control for level of sparsity before computing SHD
  G_est = RCI_res[[i]]$skel$G
  for (k in seq_len(p)){
    G_est[K[k],K[seq_len(k-1)]]=FALSE
  }
  G_est = screen_anc_Y(suffStat, G_est, SNPs=NULL, alpha=0.05)
  # # colnames(G_est) = gsub("\\..*","",colnames(ge_fs))
  CausalCell_res[[i]]$time_order = (proc.time() - ptm)[3]
  CausalCell_res[[i]]$G_est = G_est
  
  CausalCell_res[[i]]$SHD = sum(G_est!=G$graph)
  CausalCell_res[[i]]$MCC_graph = graph_MCC(G$graph>0,G_est)
 
  TWRCI_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,aa$SNPs,X,Y,B,Xte,Yte,Bte)
  locs_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,ng_1,X,Y,B,Xte,Yte,Bte)
  cis_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,cis,X,Y,B,Xte,Yte,Bte)
  ctwas_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,cTWAS,X,Y,B,Xte,Yte,Bte)
  cis_eQTLs_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,cis_eQTLs,X,Y,B,Xte,Yte,Bte)
  coloc_ABF_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,CL,X,Y,B,Xte,Yte,Bte)
  coloc_susie_res[[i]]$MACR_annot_graph_CC = MACR_annot_graph(G_est,CL_susie,X,Y,B,Xte,Yte,Bte)
  
  TWRCI_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  locs_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  cis_eQTLs_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_ABF_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  coloc_susie_res[[i]]$MSE_CRCE_CC = MSE_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,truth)
  
  TWRCI_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,aa$SNPs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  locs_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,ng_1,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,cis,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  cis_eQTLs_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,cis_eQTLs,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_ABF_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,CL,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  coloc_susie_res[[i]]$MACR_CRCE_CC = MACR_CRCE(G_est,CL_susie,X,Yg,B,Yy,Xte,Ygte,Bte,Yyte)
  
  save(file="Results_synth.RData", Gs, 
       TWRCI_res, locs_res, cis_res, ctwas_res, cis_eQTLs_res, coloc_ABF_res, coloc_susie_res, 
       SIGNET_res, RCI_res, GRCI_res, CausalCell_res)
}


# annotation
imax = 60
res_mat = array(0,c(imax,7,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MACR_annot_graph
  res_mat[i,2,1] = locs_res[[i]]$MACR_annot_graph
  res_mat[i,3,1] = cis_res[[i]]$MACR_annot_graph
  res_mat[i,4,1] = cis_eQTLs_res[[i]]$MACR_annot_graph
  res_mat[i,5,1] = coloc_ABF_res[[i]]$MACR_annot_graph
  res_mat[i,6,1] = coloc_susie_res[[i]]$MACR_annot_graph
  res_mat[i,7,1] = ctwas_res[[i]]$MACR_annot_graph      
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,4,2] = cis_eQTLs_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,5,2] = coloc_ABF_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,6,2] = coloc_susie_res[[i]]$MACR_annot_graph_SIGNET
  res_mat[i,7,2] = ctwas_res[[i]]$MACR_annot_graph_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,2,3] = locs_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,3,3] = cis_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,4,3] = cis_eQTLs_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,5,3] = coloc_ABF_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,6,3] = coloc_susie_res[[i]]$MACR_annot_graph_RCI
  res_mat[i,7,3] = ctwas_res[[i]]$MACR_annot_graph_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,4,4] = cis_eQTLs_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,5,4] = coloc_ABF_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,6,4] = coloc_susie_res[[i]]$MACR_annot_graph_GRCI
  res_mat[i,7,4] = ctwas_res[[i]]$MACR_annot_graph_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MACR_annot_graph_CC
  res_mat[i,2,5] = locs_res[[i]]$MACR_annot_graph_CC
  res_mat[i,3,5] = cis_res[[i]]$MACR_annot_graph_CC
  res_mat[i,4,5] = cis_eQTLs_res[[i]]$MACR_annot_graph_CC
  res_mat[i,5,5] = coloc_ABF_res[[i]]$MACR_annot_graph_CC
  res_mat[i,6,5] = coloc_susie_res[[i]]$MACR_annot_graph_CC
  res_mat[i,7,5] = ctwas_res[[i]]$MACR_annot_graph_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)
t(apply(res_mat,c(2,3),mean))-1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)

# CRCE
imax = 60
res_mat = array(0,c(imax,6,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MACR_CRCE
  res_mat[i,2,1] = locs_res[[i]]$MACR_CRCE
  res_mat[i,3,1] = cis_res[[i]]$MACR_CRCE
  res_mat[i,4,1] = cis_eQTLs_res[[i]]$MACR_CRCE
  res_mat[i,5,1] = coloc_ABF_res[[i]]$MACR_CRCE
  res_mat[i,6,1] = coloc_susie_res[[i]]$MACR_CRCE
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,4,2] = cis_eQTLs_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,5,2] = coloc_ABF_res[[i]]$MACR_CRCE_SIGNET
  res_mat[i,6,2] = coloc_susie_res[[i]]$MACR_CRCE_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MACR_CRCE_RCI
  res_mat[i,2,3] = locs_res[[i]]$MACR_CRCE_RCI
  res_mat[i,3,3] = cis_res[[i]]$MACR_CRCE_RCI
  res_mat[i,4,3] = cis_eQTLs_res[[i]]$MACR_CRCE_RCI
  res_mat[i,5,3] = coloc_ABF_res[[i]]$MACR_CRCE_RCI
  res_mat[i,6,3] = coloc_susie_res[[i]]$MACR_CRCE_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,4,4] = cis_eQTLs_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,5,4] = coloc_ABF_res[[i]]$MACR_CRCE_GRCI
  res_mat[i,6,4] = coloc_susie_res[[i]]$MACR_CRCE_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MACR_CRCE_CC
  res_mat[i,2,5] = locs_res[[i]]$MACR_CRCE_CC
  res_mat[i,3,5] = cis_res[[i]]$MACR_CRCE_CC
  res_mat[i,4,5] = cis_eQTLs_res[[i]]$MACR_CRCE_CC
  res_mat[i,5,5] = coloc_ABF_res[[i]]$MACR_CRCE_CC
  res_mat[i,6,5] = coloc_susie_res[[i]]$MACR_CRCE_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)
t(apply(res_mat,c(2,3),mean))-1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)


# MSE CRCE
imax = 60
res_mat = array(0,c(imax,6,5))
for (i in 1:imax){
  res_mat[i,1,1] = TWRCI_res[[i]]$MSE_CRCE
  res_mat[i,2,1] = locs_res[[i]]$MSE_CRCE
  res_mat[i,3,1] = cis_res[[i]]$MSE_CRCE
  res_mat[i,4,1] = cis_eQTLs_res[[i]]$MSE_CRCE
  res_mat[i,5,1] = coloc_ABF_res[[i]]$MSE_CRCE
  res_mat[i,6,1] = coloc_susie_res[[i]]$MSE_CRCE
  
  res_mat[i,1,2] = TWRCI_res[[i]]$MSE_CRCE_SIGNET
  res_mat[i,2,2] = locs_res[[i]]$MSE_CRCE_SIGNET
  res_mat[i,3,2] = cis_res[[i]]$MSE_CRCE_SIGNET
  res_mat[i,4,2] = cis_eQTLs_res[[i]]$MSE_CRCE_SIGNET
  res_mat[i,5,2] = coloc_ABF_res[[i]]$MSE_CRCE_SIGNET
  res_mat[i,6,2] = coloc_susie_res[[i]]$MSE_CRCE_SIGNET
  
  res_mat[i,1,3] = TWRCI_res[[i]]$MSE_CRCE_RCI
  res_mat[i,2,3] = locs_res[[i]]$MSE_CRCE_RCI
  res_mat[i,3,3] = cis_res[[i]]$MSE_CRCE_RCI
  res_mat[i,4,3] = cis_eQTLs_res[[i]]$MSE_CRCE_RCI
  res_mat[i,5,3] = coloc_ABF_res[[i]]$MSE_CRCE_RCI
  res_mat[i,6,3] = coloc_susie_res[[i]]$MSE_CRCE_RCI
  
  res_mat[i,1,4] = TWRCI_res[[i]]$MSE_CRCE_GRCI
  res_mat[i,2,4] = locs_res[[i]]$MSE_CRCE_GRCI
  res_mat[i,3,4] = cis_res[[i]]$MSE_CRCE_GRCI
  res_mat[i,4,4] = cis_eQTLs_res[[i]]$MSE_CRCE_GRCI
  res_mat[i,5,4] = coloc_ABF_res[[i]]$MSE_CRCE_GRCI
  res_mat[i,6,4] = coloc_susie_res[[i]]$MSE_CRCE_GRCI
  
  res_mat[i,1,5] = TWRCI_res[[i]]$MSE_CRCE_CC
  res_mat[i,2,5] = locs_res[[i]]$MSE_CRCE_CC
  res_mat[i,3,5] = cis_res[[i]]$MSE_CRCE_CC
  res_mat[i,4,5] = cis_eQTLs_res[[i]]$MSE_CRCE_CC
  res_mat[i,5,5] = coloc_ABF_res[[i]]$MSE_CRCE_CC
  res_mat[i,6,5] = coloc_susie_res[[i]]$MSE_CRCE_CC
  
}

t(apply(res_mat,c(2,3),mean))
t(apply(res_mat,c(2,3),mean))+1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)
t(apply(res_mat,c(2,3),mean))-1.96*t(apply(res_mat,c(2,3),sd))/sqrt(60)


# annotation
imax = 60
res_mat = matrix(0,imax,7)
for (i in 1:imax){
  res_mat[i,1] = TWRCI_res[[i]]$MCC
  res_mat[i,2] = locs_res[[i]]$MCC
  res_mat[i,3] = cis_res[[i]]$MCC
  res_mat[i,4] = cis_eQTLs_res[[i]]$MCC
  res_mat[i,5] = coloc_ABF_res[[i]]$MCC
  res_mat[i,6] = coloc_susie_res[[i]]$MCC
  res_mat[i,7] = ctwas_res[[i]]$MCC
}
print(colMeans(res_mat))

# graph
imax = 60
res_mat = matrix(0,imax,5)
for (i in 1:imax){
  res_mat[i,1] = TWRCI_res[[i]]$SHD
  res_mat[i,2] = SIGNET_res[[i]]$SHD
  res_mat[i,3] = RCI_res[[i]]$SHD
  res_mat[i,4] = GRCI_res[[i]]$SHD
  res_mat[i,5] = CausalCell_res[[i]]$SHD
}
print(colMeans(res_mat))

# timing annotation
imax = 60
res_mat = matrix(0,imax,7)
for (i in 1:imax){
  print(i)
  res_mat[i,1] = TWRCI_res[[i]]$time
  res_mat[i,2] = 0
  res_mat[i,3] = cis_res[[i]]$time
  res_mat[i,4] = cis_eQTLs_res[[i]]$time
  res_mat[i,5] = coloc_ABF_res[[i]]$time
  res_mat[i,6] = coloc_susie_res[[i]]$time
  res_mat[i,7] = ctwas_res[[i]]$time
}
print(colMeans(res_mat))

# timing graph
imax = 60
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

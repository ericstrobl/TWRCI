ctwas2 <- function(Z,X,Y,cis,batches){

  require(logging)
  require(ctwas)
  require(earth)
  
  X = earth(batches,X)$residuals
  Y = earth(batches,Y)$residuals
  
  Z = normalizeData(Z)
  X = normalizeData(X)
  Y = normalizeData(Y)
  
  n = nrow(Z)
  p = ncol(Z)
  
  ## SNP stuff, put all on chr 16
  r_snp = cor(Z,Y)
  z_snp = r_snp*sqrt((n-2)/(1-r_snp^2))
  colnames(z_snp) = c()
  z_snp = data.frame(id = 1:p,A1 = rep("A",p),A2 = rep("C",p),z=z_snp)

  R_snp = cor(Z)
  saveRDS(R_snp,file="~/summary_graph/ctwas_stuff/R_snp.RDS")
  
  R_snp_info = data.frame(chrom = rep(16,p),id=1:p,pos=1:p,alt=rep("C",p),ref=rep("A",p),
                    variance = apply(Z,2,var),allele_freq = apply(Z>0,2,mean))
  write.table(R_snp_info,file="~/summary_graph/ctwas_stuff/R_snp.Rvar",sep = '\t',row.names=F,quote=F)
  
  ld_R_dir = "~/summary_graph/ctwas_stuff"
  
  regions = data.frame(chrom = 16,start=1,stop=(p+1))
  write.table(regions,file="~/summary_graph/ctwas_stuff/regions_subset.bed",sep = '\t',row.names=F,quote=F)
  regions_file = "~/summary_graph/ctwas_stuff/regions_subset.bed"
  
  ## imputation stuff
  X_imp = c()
  d = ncol(X)
  wgtlist = list()
  ig = c(); j =1
  for (i in 1:d){
    if(length(cis[[i]])>0){
      mod = CV_LRR2_all(Z[,cis[[i]]],X[,i])
      X_imp = cbind(X_imp,mod$fitted.values)
      rownames(mod$beta0) = cis[[i]]
      colnames(mod$beta0) = "weight"
      wgtlist[[j]] = mod$beta0
      j = j+1
      ig = c(ig, i)
    } #else{
      #X_imp[,i] = mean(X)
      #wgtlist[i] = list(NULL)
    #}
  }
  names(wgtlist)=paste0("G",ig,sep="")
  save(wgtlist,file="~/summary_graph/ctwas_stuff/example_locus_chr16.exprqc.Rd")
  
  # is = which(apply(X_imp,2,sd)>0)
  r_gene = cor(X_imp,Y)
  z_gene = r_gene*sqrt((n-2)/(1-r_gene^2))
  colnames(z_gene)=c()
  z_gene = data.frame(id = paste("G",ig,sep=""),z=z_gene)
  
  pos = t(matrix(1:(2*length(ig)),nrow=2))
  chr16 = data.frame(chrom = rep(16,length(ig)),id=paste("G",ig,sep=""),p0=pos[,1],p1=pos[,2])
  write.table(chr16,file="C:/users/ericv/Documents/summary_graph/ctwas_stuff/example_locus_chr16.exprvar",sep = '\t',row.names=F,quote=F)
  
  imp_LD = data.frame(chrom = rep(16,length(ig)),start=pos[,1],stop=pos[,2],
                      RDS_file= "C:/users/ericv/Documents/summary_graph/ctwas_stuff/R_snp.RDS" ,
                        region_name = rep(1,length(ig)))
  write.table(imp_LD,file="~/summary_graph/ctwas_stuff/example_locus_ld_R_chr16.txt",sep = '\t',row.names=F,quote=F)
  
  outputdir = "~/summary_graph/ctwas_stuff/output"
  outname = "CTWAS"
  
  ld_exprvarfs =  paste0("~/summary_graph/ctwas_stuff/example_locus_chr", 1:22, ".exprvar")
  
  if (sum(!is.na(z_gene$z))==0){
    z_gene = NULL
  }
  ctwas_rss2(z_gene = z_gene,
             z_snp = z_snp,
             ld_exprvarfs = ld_exprvarfs,
             ld_R_dir = ld_R_dir,
             ld_regions_custom = regions_file,
             outputdir = outputdir,
             outname = outname)
  
  ctwas_res <- read.table(paste0(outputdir, "/", outname, ".susieIrss.txt"), header=T)
  iS = which(ctwas_res$type == "SNP")
  cis[[d+1]] = which( ctwas_res$susie_pip[iS]>0.8 )

  return(cis)
}


ctwas_rss2 <- function (z_gene, z_snp, ld_exprvarfs, ld_exprfs = NULL, ld_pgenfs = NULL, 
          ld_R_dir = NULL, ld_regions = c("EUR", "ASN", "AFR"), ld_regions_version = c("b37", 
                                                                                       "b38"), ld_regions_custom = NULL, thin = 1, prob_single = 0.0, 
          rerun_gene_PIP = 0.8, niter1 = 3, niter2 = 30, L = 5, group_prior = NULL, 
          group_prior_var = NULL, estimate_group_prior = FALSE, estimate_group_prior_var = FALSE, 
          use_null_weight = TRUE, coverage = 0.95, max_snp_region = Inf, 
          ncore = 1, ncore.rerun = 1, outputdir = getwd(), outname = NULL, 
          logfile = NULL, merge = TRUE) 
{
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  logging:::loginfo("ctwas started ... ")
  if (length(ld_exprvarfs) != 22) {
    stop("Not all imputed expression files for 22 chromosomes are provided.")
  }
  if (is.null(ld_pgenfs) & is.null(ld_R_dir)) {
    stop("Error: need to provide either .pgen file or ld_R file")
  }
  else if (!is.null(ld_pgenfs)) {
    if (length(ld_pgenfs) != 22) {
      stop("Not all pgen files for 22 chromosomes are provided.")
    }
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_pvarfs, ctwas:::read_pvar)
    ld_Rfs <- NULL
  }
  else {
    ld_Rfs <- ctwas:::write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_Rfs, ctwas:::read_ld_Rvar)
    ld_pvarfs <- NULL
  }
  if (is.null(ld_regions_custom)) {
    ld_regions <- match.arg(ld_regions)
    ld_regions_version <- match.arg(ld_regions_version)
    regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, 
                                                           ".", ld_regions_version, ".bed"), package = "ctwas")
  }
  else {
    regionfile <- ld_regions_custom
  }
  logging:::loginfo("LD region file: %s", regionfile)
  zdf <- rbind(z_snp[, c("id", "z")], z_gene[, c("id", "z")])
  rm(z_snp, ld_snpinfo)
  if (thin <= 0 | thin > 1) {
    stop("thin value needs to be in (0,1]")
  }
  regionlist <- ctwas:::index_regions(regionfile = regionfile, exprvarfs = ld_exprvarfs, 
                              pvarfs = ld_pvarfs, ld_Rfs = ld_Rfs, select = zdf$id, 
                              thin = thin, minvar = 2, outname = outname, outputdir = outputdir, 
                              merge = merge)
  saveRDS(regionlist, file = paste0(outputdir, "/", outname, 
                                    ".regionlist.RDS"))
  temp_regs <- lapply(1:22, function(x) cbind(x, unlist(lapply(regionlist[[x]], 
                                                               "[[", "start")), unlist(lapply(regionlist[[x]], "[[", 
                                                                                              "stop"))))
  regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 
                                                           3) {
    x
  }))
  write.table(regs, file = paste0(outputdir, "/", outname, 
                                  ".regions.txt"), row.names = F, sep = "\t", 
              quote = F)
  if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)) {
    logging:::loginfo("Run susie iteratively, getting rough estimate ...")
    if (!is.null(group_prior)) {
      group_prior[2] <- group_prior[2]/thin
    }

    pars <- susieI_rss2(zdf = zdf, regionlist = regionlist, 
                       ld_exprvarfs = ld_exprvarfs, ld_exprfs = ld_exprfs, 
                       ld_pgenfs = ld_pgenfs, ld_Rfs = ld_Rfs, niter = niter1, 
                       L = 1, z_ld_weight = 0, group_prior = group_prior, 
                       group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                       estimate_group_prior_var = estimate_group_prior_var, 
                       use_null_weight = use_null_weight, coverage = coverage, 
                       ncore = ncore, outputdir = outputdir, outname = paste0(outname, 
                                                                              ".s1"))
    group_prior <- pars[["group_prior"]]
    group_prior_var <- pars[["group_prior_var"]]
    regionlist2 <- ctwas:::filter_regions(regionlist, group_prior, 
                                  prob_single = prob_single)
    logging:::loginfo("Blocks are filtered: %s blocks left", sum(unlist(lapply(regionlist2, 
                                                                     length))))
    logging:::loginfo("Run susie iteratively, getting accurate estimate ...")
    pars <- susieI_rss2(zdf = zdf, regionlist = regionlist2, 
                       ld_exprvarfs = ld_exprvarfs, ld_exprfs = ld_exprfs, 
                       ld_pgenfs = ld_pgenfs, ld_Rfs = ld_Rfs, niter = niter2, 
                       L = 1, z_ld_weight = 0, group_prior = group_prior, 
                       group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                       estimate_group_prior_var = estimate_group_prior_var, 
                       use_null_weight = use_null_weight, coverage = coverage, 
                       ncore = ncore, outputdir = outputdir, outname = paste0(outname, 
                                                                              ".s2"))
    group_prior <- pars[["group_prior"]]
    group_prior_var <- pars[["group_prior_var"]]
  }
  logging:::loginfo("Run susie for all regions.")
  pars <- susieI_rss2(zdf = zdf, regionlist = regionlist, ld_exprvarfs = ld_exprvarfs, 
                     ld_exprfs = ld_exprfs, ld_pgenfs = ld_pgenfs, ld_Rfs = ld_Rfs, 
                     niter = 1, L = L, z_ld_weight = 0, group_prior = group_prior, 
                     group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                     estimate_group_prior_var = estimate_group_prior_var, 
                     use_null_weight = use_null_weight, coverage = coverage, 
                     ncore = ncore, outputdir = outputdir, outname = paste0(outname, 
                                                                            ".temp"), report_parameters = F)
  group_prior[2] <- group_prior[2] * thin
  if (thin == 1) {
    file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"), 
                paste0(file.path(outputdir, outname), ".susieIrss.txt"))
  }
  else {
    regionlist <- index_regions(regionfile = regionfile, 
                                exprvarfs = ld_exprvarfs, pvarfs = ld_pvarfs, ld_Rfs = ld_Rfs, 
                                select = zdf, thin = 1, maxSNP = max_snp_region, 
                                minvar = 2, outname = outname, outputdir = outputdir, 
                                merge = merge)
    res <- data.table::fread(paste0(file.path(outputdir, 
                                              outname), ".temp.susieIrss.txt"))
    res.keep <- NULL
    for (b in 1:length(regionlist)) {
      for (rn in names(regionlist[[b]])) {
        gene_PIP <- max(res[res$type == "gene" & res$region_tag1 == 
                              b & res$region_tag2 == rn, ]$susie_pip, 0)
        if (gene_PIP < rerun_gene_PIP) {
          regionlist[[b]][[rn]] <- NULL
          res.keep <- rbind(res.keep, res[res$region_tag1 == 
                                            b & res$region_tag2 == rn, ])
        }
      }
    }
    nreg <- sum(unlist(lapply(regionlist, length)))
    logging:::loginfo("Number of regions that contain strong gene signals: %s", 
            nreg)
    if (nreg == 0) {
      file.rename(paste0(file.path(outputdir, outname), 
                         ".temp.susieIrss.txt"), paste0(file.path(outputdir, 
                                                                  outname), ".susieIrss.txt"))
    }
    else {
      logging:::loginfo("Rerun susie for regions with strong gene signals using full SNPs.")
      pars <- susieI_rss2(zdf = zdf, regionlist = regionlist, 
                         ld_exprvarfs = ld_exprvarfs, ld_exprfs = ld_exprfs, 
                         ld_pgenfs = ld_pgenfs, ld_Rfs = ld_Rfs, niter = 1, 
                         L = L, z_ld_weight = 0, group_prior = group_prior, 
                         group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                         estimate_group_prior_var = estimate_group_prior_var, 
                         use_null_weight = use_null_weight, coverage = coverage, 
                         ncore = ncore.rerun, outputdir = outputdir, outname = paste0(outname, 
                                                                                      ".s3"), report_parameters = F)
      res.rerun <- data.table::fread(paste0(file.path(outputdir, 
                                                      outname), ".s3.susieIrss.txt"))
      res <- rbind(res.keep, res.rerun)
      data.table::fwrite(res, file = paste0(file.path(outputdir, 
                                                      outname), ".susieIrss.txt"), sep = "\t", quote = F)
      file.remove((paste0(file.path(outputdir, outname), 
                          ".temp.susieIrss.txt")))
    }
  }
  list(group_prior = group_prior, group_prior_var = group_prior_var)
}

anno_susie2 <- function (susieres, geneinfo, snpinfo, gidx, sidx, region_tag1, 
          region_tag2) 
{
  anno.gene <- NULL
  if (length(geneinfo) != 0) {
    anno.gene <- cbind(geneinfo[gidx, c("chrom", "id", "p0")], 
                       rep("gene", length(gidx)))
    colnames(anno.gene) <- c("chrom", "id", "pos", "type")
  }
  anno.SNP <- cbind(snpinfo[sidx, c("chrom", "id", "pos")], 
                    rep("SNP", length(sidx)))
  colnames(anno.SNP) <- c("chrom", "id", "pos", "type")
  anno <- rbind(anno.gene, anno.SNP)
  anno <- as.data.frame(anno)
  anno$region_tag1 <- region_tag1
  anno$region_tag2 <- region_tag2
  anno$cs_index <- 0
  if (!is.null(susieres$sets$cs)) {
    for (cs_i in susieres$sets$cs_index) {
      X.idx <- susieres$sets$cs[[paste0("L", cs_i)]]
      X.idx <- X.idx[X.idx != susieres$null_index]
      anno$cs_index[X.idx] <- cs_i
    }
  }
  outdf.rn <- cbind(anno, susieres$pip)
  colnames(outdf.rn)[8] <- "susie_pip"
  p <- length(gidx) + length(sidx)
  outdf.rn$mu <- colSums(susieres$mu[, seq(1, p)[1:p != susieres$null_index], 
                                       drop = F])
  outdf.rn$mu2 <- colSums(susieres$mu2[, seq(1, p)[1:p != susieres$null_index], 
                                       drop = F])
  outdf.rn
}

susieI_rss2 <- function (zdf, regionlist, ld_exprvarfs, ld_exprfs = NULL, ld_pgenfs = NULL, 
          ld_Rfs = NULL, niter = 20, L = 1, z_ld_weight = 0, group_prior = NULL, 
          group_prior_var = NULL, estimate_group_prior = TRUE, estimate_group_prior_var = TRUE, 
          use_null_weight = TRUE, coverage = 0.95, ncore = 1, outputdir = getwd(), 
          outname = NULL, report_parameters = TRUE) 
{
  require(foreach)
  outname <- file.path(outputdir, outname)
  if (is.null(ld_pgenfs) & is.null(ld_Rfs)) {
    stop("Error: need to provide either .pgen file or ld_R file")
  }
  if (!is.null(ld_pgenfs)) {
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
  }
  K <- 2
  group_prior_rec <- matrix(, nrow = K, ncol = niter)
  group_prior_var_rec <- matrix(, nrow = K, ncol = niter)
  prior.gene <- group_prior[1]
  prior.SNP <- group_prior[2]
  V.gene <- group_prior_var[1]
  V.SNP <- group_prior_var[2]
  for (iter in 1:niter) {
    logging:::loginfo("run iteration %s", iter)
    snp.rpiplist <- list()
    gene.rpiplist <- list()
    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)
    corelist <- ctwas:::region2core(regionlist, ncore)
    outdf <- foreach(core = 1:length(corelist), .combine = "rbind", .export = c("anno_susie2"),
                     .packages = "ctwas") %dopar% {
                       outdf.core.list <- list()
                       regs <- corelist[[core]]
                       for (reg in 1:nrow(regs)) {
                         b <- regs[reg, "b"]
                         rn <- regs[reg, "rn"]
                         gidx <- regionlist[[b]][[rn]][["gidx"]]
                         sidx <- regionlist[[b]][[rn]][["sidx"]]
                         gid <- regionlist[[b]][[rn]][["gid"]]
                         sid <- regionlist[[b]][[rn]][["sid"]]
                         p <- length(gidx) + length(sidx)
                         if (is.null(prior.gene) | is.null(prior.SNP)) {
                           prior <- c(rep(1/p, length(gidx)), rep(1/p, 
                                                                  length(sidx)))
                         }
                         else {
                           prior <- c(rep(prior.gene, length(gidx)), rep(prior.SNP, 
                                                                         length(sidx)))
                         }
                         if (is.null(V.gene) | is.null(V.SNP)) {
                           V <- matrix(rep(50, L * p), nrow = L)
                         }
                         else {
                           V <- c(rep(V.gene, length(gidx)), rep(V.SNP, 
                                                                 length(sidx)))
                           V <- matrix(rep(V, each = L), nrow = L)
                         }
                         if (isTRUE(use_null_weight)) {
                           nw <- max(0, 1 - sum(prior))
                           prior <- prior/(1 - nw)
                         }
                         else {
                           nw <- NULL
                         }
                         z.g <- zdf[match(gid, zdf$id), ][["z"]]
                         z.s <- zdf[match(sid, zdf$id), ][["z"]]
                         z <- c(z.g, z.s)
                         if (!(is.null(ld_pgenfs))) {
                           ld_pgen <- ctwas:::prep_pgen(pgenf = ld_pgenfs[b], 
                                                ld_pvarfs[b])
                           X.g <- ctwas:::read_expr(ld_exprfs[b], variantidx = gidx)
                           X.s <- ctwas:::read_pgen(ld_pgen, variantidx = sidx)
                           X <- cbind(X.g, X.s)
                           R <- Rfast::cora(X)
                         }
                         else {
                           if (!("regRDS" %in% names(regionlist[[b]][[rn]]))) {
                             stop("R matrix info not available for region", 
                                  b, ",", rn)
                           }
                           regRDS <- regionlist[[b]][[rn]][["regRDS"]]
                           if (is.null(regionlist[[b]][[rn]][["R_s_file"]])) {
                             R_snp <- lapply(regRDS, readRDS)
                             R_snp <- suppressWarnings({
                               as.matrix(Matrix::bdiag(R_snp))
                             })
                             R_snp <- R_snp[sidx, sidx, drop = F]
                           }
                           else {
                             R_snp <- readRDS(regionlist[[b]][[rn]][["R_s_file"]])
                           }
                           R_snp_gene <- readRDS(regionlist[[b]][[rn]][["R_sg_file"]])
                           R_snp_gene <- R_snp_gene[sidx, , drop = F]
                           R_gene <- readRDS(regionlist[[b]][[rn]][["R_g_file"]])
                           R <- rbind(cbind(R_gene, t(R_snp_gene)), cbind(R_snp_gene, 
                                                                          R_snp))
                         }
 
                         susieres <- ctwas:::susie_rss(z, R, z_ld_weight = z_ld_weight, 
                                               L = L, prior_weights = prior, null_weight = nw, 
                                               prior_variance = V, estimate_prior_variance = F, 
                                               coverage = coverage)
                         geneinfo <- ctwas:::read_exprvar(ld_exprvarfs[b])
                         if (!is.null(ld_pgenfs)) {
                           snpinfo <- ctwas:::read_pvar(ld_pvarfs[b])
                         }
                         else {
                           snpinfo <- do.call(rbind, lapply(regRDS, ctwas:::read_ld_Rvar_RDS))
                         }
                         outdf.reg <- anno_susie2(susieres, geneinfo, snpinfo, 
                                                 gidx, sidx, b, rn)
                         outdf.core.list[[reg]] <- outdf.reg
                       }
                       outdf.core <- do.call(rbind, outdf.core.list)
                       outdf.core
                     }
    if (isTRUE(estimate_group_prior)) {
      prior.SNP <- mean(outdf[outdf[, "type"] == "SNP", 
                              "susie_pip"])
      prior.gene <- mean(outdf[outdf[, "type"] == "gene", 
                               "susie_pip"])
      group_prior_rec[, iter] <- c(prior.gene, prior.SNP)
    }
    if (report_parameters) {
      logging:::loginfo("After iteration %s, gene prior %s:, SNP prior:%s", 
              iter, prior.gene, prior.SNP)
    }
    if (isTRUE(estimate_group_prior_var)) {
      outdf.g <- outdf[outdf[, "type"] == "gene", ]
      outdf.s <- outdf[outdf[, "type"] == "SNP", ]
      V.gene <- sum(outdf.g$susie_pip * outdf.g$mu2)/sum(outdf.g$susie_pip)
      V.SNP <- sum(outdf.s$susie_pip * outdf.s$mu2)/sum(outdf.s$susie_pip)
      group_prior_var_rec[, iter] <- c(V.gene, V.SNP)
    }
    save(group_prior_rec, group_prior_var_rec, file = paste0(outname, 
                                                             ".susieIrssres.Rd"))
    data.table::fwrite(outdf, file = paste0(outname, ".susieIrss.txt"), 
                       sep = "\t", quote = F)
    parallel::stopCluster(cl)
  }
  list(group_prior = c(prior.gene, prior.SNP), group_prior_var = c(V.gene, 
                                                                   V.SNP))
}

                                

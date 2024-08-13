# Transcriptome-Wide Root Causal Inference (TWRCI)

This is an R package implementing the TWRCI algorithm for discovering root causal gene expression levels -- or root causal genes for short -- from observational data. Root causal genes correspond to the first gene expression levels that are disturbed in disease near the beginning of pathogenesis. In contrast, core genes lie at the end of pathogenesis and driver genes only account for the effects of somatic mutations primarily in cancer. TWRCI requires individual level data containing genetic variants, gene expression levels from the relevant tissue and the phenotype (variant-expression-phenotype data).

The academic article describing TWRCI in detail can be found [here](https://www.medrxiv.org/content/10.1101/2024.07.22.24310837v2). Please cite the article if you use any of the code in this repository.

The Experiments folder contains any additional code needed to replicate the experimental results in the paper after downloading [GTEx V8 protected access data](https://gtexportal.org/home/protectedDataAccess) and lifting over to hg19 with [BCFtools](https://samtools.github.io/bcftools/). All code was tested in R version 4.3.1 and BCFtools version 1.18.

# Installation
Install time is under 5 minutes for a typical desktop computer.

> if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> BiocManager::install(c("qvalue", "biomaRt", "AnnotationDbi", "org.Hs.eg.db", "gage"))

> install_github("chaim-e/independence")

> library(devtools)

> install_github("ericstrobl/TWRCI")

> library(TWRCI)

Installation of the cTWAS package is also required if you want to replicate the experiments:

> install_github("xinhe-lab/ctwas",ref = "main")

# Run TWRCI on demo data
Set the number of samples of the variant-expression-phenotype data to 500:

> nsamps = 500

Instantiate 30 instrumental variables:

> datat = matrix(rnorm(nsamps*30),nsamps,30)

Set the number of gene expression levels plus the phenotype to 10:

> p = 10 

Assign variants to each gene and the phenotype:

> SNPs = generate_SNP_list(p,ncol(datat))

Generate the DAG:

> G = generate_DAG(p,SNPs$SNPs,bulk_nbatch=1)

Generate the remaining expression-phenotype data:

> RNAseq = sample_data(G,SNPs$SNPs,datat,p=p,nsamps=nsamps)

Run TWRCI:

> out = TWRCI(RNAseq$X, RNAseq$SNP_data, RNAseq$Y, RNAseq$batch, RNAseq$batch, graph=TRUE, CRCEs=TRUE)

Run time is under 2 minutes for a typical desktop computer. Print the output of TWRCI:

> print(out$K) # causal order

> print(out$SNPs) # annotations

> print(out$G_est) # causal graph

> print(out$CRCEs) # CRCE estimates 


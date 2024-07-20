# Transcriptome-Wide Root Causal Inference (TWRCI)

This is an R package implementing the TWRCI algorithm for discovering root causal gene expression levels -- or root causal genes for short -- from observational data. Root causal genes correspond to the first gene expression levels that are disturbed in disease and help initiate pathogenesis. In contrast, core genes lie at the end of pathogenesis and driver genes only account for the effects of somatic mutations primarily in cancer. TWRCI requires individual level data containing genetic variants, gene expression levels from the relevant tissue and the phenotype (variant-expression-phenotype data).

The academic article describing TWRCI in detail can be found [here](https://www.google.com). Please cite the article if you use any of the code in this repository.

The Experiments folder contains any additional code needed to replicate the experimental results in the paper. All code was tested in R version 4.3.1.

# Installation
> if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> library(devtools)

> install_github("ericstrobl/TWRCI")

> library(TWRCI)

# Run TWRCI on sample data
Set the number of samples of the variant-expression-phenotype data to 500:

> nsamps = 500

Instantiate 100 instrumental variables:

> datat = matrix(rnorm(nsamps*100),nsamps,100)

Set the number of gene expression levels plus the phenotype to 30:

> p = 30 

Assign SNPs to each gene:

> SNPs = generate_SNP_list(p,ncol(data)/p)

Generate the DAG:

> G = generate_DAG(p,SNPs$SNPs,bulk_nbatch=1)

Generate the reamining expression-phenotype data:

> RNAseq = sample_data(G,SNPs$SNPs,datat,p=p,nsamps=nsamps)

Run TWRCI:

> out = TWRCI(RNAseq$X,RNAseq$SNP_data,RNAseq$Y,RNAseq$batch)

Print output of TWRCI:
> print(out$K) # causal order
> print(out$SNPs) # annotations
> print(out$G_est) # causal graph
> print(out$CRCEs) # CRCE estimates 


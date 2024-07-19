# Transcriptome-Wide Root Causal Inference (TWRCI)

This is an R package implementing the TWRCI algorithm for discovering root causal gene expression levels -- or root causal genes for short -- from observational data. Root causal genes correspond to the first gene expression levels that are disturbed in disease and help initiate pathogenesis. In contrast, core genes lie at the end of pathogenesis and driver genes only account for the effects of somatic mutations primarily in cancer. TWRCI requires individual level data containing genetic variants, gene expression levels from the relevant tissue and the phenotype (variant-expression-phenotype data).

The academic article describing TWRCI in detail can be found [here](https://www.google.com). Please cite the article if you use any of the code in this repository.

The Experiments folder contains any additional code needed to replicate the experimental results in the paper. All code was tested in R version 4.3.1.

# Installation
> if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> BiocManager::install(c("fgsea","qvalue","org.Hs.eg.db","reactome.db"))

> library(devtools)

> install_github("ericstrobl/TWRCI")

> library(TWRCI)

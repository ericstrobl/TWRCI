### run this script after lifting the GTEx individual-level variant data with BCFtools to hg19 (tested on version 1.18)

require(ieugwasr)

leads = tophits("finn-b-I9_ISCHHEART",pval = 5e-05, clump = 0)

ts = abs(qt(leads$p/2,leads$n-2))
leads = cbind(leads,r=ts/sqrt(leads$n+ts^2-2)*sign(leads$beta)) # compute correlation coefficient

###########
# gene expression data
ge = read.csv("GTEx_gene_reads_whole_blood.csv",header=TRUE)
colnames = ge[,2]
ge = ge[,-(1:3)]
ge = t(ge)
colnames(ge) = colnames
ix = which(colMeans(ge)<=5)
ge = ge[,-ix,drop=FALSE]
rownames(ge) = sub("^(([^.]*.){1}[^.]*).*", "\\1", rownames(ge))
save(file="GTEx_gene_reads_whole_blood_preprocessed.RData",ge)

#####
# SNP data
library(gwasvcf)
chrpos = paste(leads$chr,":",leads$position,sep="")

vcf_list = my_query_chrompos_file(chrpos, vcffile='C:/cygwin64/home/ericv/bin/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.sorted.vcf.gz')
data = vcf_list$vcf
chrompos = vcf_list$chrompos
SNPn = gt2num(data@assays@data$GT)
SNP_data = t(SNPn$genomat)
rownames(SNP_data) = gsub("-",".",data@metadata$header@samples)
io=findOverlaps(data,chrompos,select="first")
colnames(SNP_data) = leads$rsid[io]

# ensure REF and ALT alleles match exactly
ielim = union(which(as.character(data@fixed$REF)!=leads$nea[io]),which(as.character(unlist(data@fixed$ALT))!=leads$ea[io])) # match REF and ALT alleles
ielim = union(ielim, which(is.na(colSums(SNP_data)))) #no NA values
ielim = union(ielim,which(duplicated(io))) #not duplicated
SNP_data = SNP_data[,-ielim]
leads_f = leads[io,][-ielim,]

save(file="IHD_whole_blood_SNP_data_preprocessed.RData",SNP_data,leads_f)


### intersect data
# samps
samps = intersect(rownames(ge),rownames(SNP_data))
SNP_data = SNP_data[which(rownames(SNP_data) %in% samps),]
ge = ge[which(rownames(ge) %in% samps),]
samps = intersect(which(complete.cases(SNP_data)),which(complete.cases(ge)))
SNP_data_f = SNP_data[samps,,drop=FALSE]
ge_f = ge[samps,,drop=FALSE]

### control variables
# phenotype
phenot = read.csv("C:/Users/ericv/Documents/summary_graph/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.csv")
phenot = t(phenot); phenot = phenot[-1,]
ID = rownames(phenot)
phenot = matrix(as.numeric(phenot),ncol=ncol(phenot))
pheno = phenot[,c(1:5,(ncol(phenot)-2):ncol(phenot))]; rownames(pheno) = ID
pheno_f = pheno[rownames(pheno) %in% rownames(ge_f),]

targett = read.csv("C:/Users/ericv/Documents/gtex data/pheno/GTEx_phenotype.csv")
target = as.matrix(targett$MHHRTDIS); #####
age = as.matrix(targett$AGE)

rownames(target) = targett[,1]
rownames(age) = targett[,1]
target = target[complete.cases(target),,drop=FALSE]
age = age[complete.cases(target),,drop=FALSE]

i_int = intersect(rownames(target),rownames(ge_f))

target_f = target[which(rownames(target) %in% i_int),,drop=FALSE]
pheno_f = pheno[which(rownames(pheno) %in% i_int),]
SNP_data_f = SNP_data_f[which(rownames(SNP_data_f) %in% i_int),]
ge_f = ge_f[which(rownames(ge_f) %in% i_int),]
age_f = age[which(rownames(age) %in% i_int),]

save(file="IHD_whole_blood_exp_SNP_data2.RData",SNP_data_f,ge_f,leads_f,pheno_f,target_f,age_f)
# pheno_f = pheno_f[,8]

ge_f = normalizeData(asinh(ge_f))
pheno_f = normalizeData(pheno_f)
SNP_data_f = normalizeData(SNP_data_f)

ge_f = ge_f - lin_reg(pheno_f,ge_f)$fitted.values
SNP_data_f = SNP_data_f - lin_reg(pheno_f,SNP_data_f)$fitted.values
target_f = target_f - lin_reg(pheno_f,target_f)$fitted.values

ge_f = normalizeData(ge_f)
SNP_data_f = normalizeData(SNP_data_f)
target_f = normalizeData(target_f)

save(file="IHD_whole_blood_exp_SNP_data_controlled2.RData",SNP_data_f,ge_f,leads_f,pheno_f,target_f,age_f)


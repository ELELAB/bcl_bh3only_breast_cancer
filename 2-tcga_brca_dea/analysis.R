#!/usr/bin/env Rscript
setwd(getwd())
start_time <- Sys.time()
source('src/functions.R')
#source('src/download_data.R') uncomment to allow of downloading data in pipeline

#### Analyses of the Cancer Genome Atlas (TCGA) RNASeq BRCA dataset  ####

# Install and load packages if not installed. Othwewise load.
source('https://bioconductor.org/biocLite.R')
list.of.packages <- c('ggplot2', 'RColorBrewer','gplots','SummarizedExperiment', 'limma', 'TCGAbiolinks',
                      'igraph','openxlsx', 'devtools', 'biomaRt', 'psych', 'dplyr', 'cgdsr')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) biocLite(new.packages)
sapply(list.of.packages, require, character.only = TRUE)
install_github('vqv/ggbiplot')
library(ggbiplot)


BCL2_genes <- c('BCL2','BCL2L1','BCL2L2','MCL1', 'BCL2L10', 'BCL2A1', 'BAX',
                'BOK','BAK1','BCL2L12',
                'BCL2L13','BCL2L14','BCL2L15')

# Candidate genes: BCl-2 homologs and interactors from the the IID mammary gland database
candidate_genes <- readLines('data/candidate.genes.txt')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### DATA #### 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load RangedSummarizedExperiment of the TCGA BRCA dataset
rse <- get(load('data/brcaExp_harmonized_2017-11-24.rda'))
#rse <- data_harmo # only used when sourcing the download script
which(duplicated(colnames(rse)))  # test of identical samples exist


clinic_df <- get(load('data/brca_clin.RData'))

# get subtype information 
data_subt <- TCGAquery_subtype(tumor = 'BRCA')

# clinical data
#clinic_df <- query_clin # only used when sourcing the download_data.R


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### PRE-PROCESSING AND NORMALIZATION #### 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Filter the data, so that it does only contain protein coding transcripts

# get list of all protein coding genes
mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
listOfGenes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 'gene_biotype'),filters = c('biotype'),values = list(biotype='protein_coding'), mart = mart)
head(listOfGenes)
rse_coding <- rse[rownames(rse) %in% listOfGenes$ensembl_gene_id, ]


# Estimate the pearson correlation coefficient among all pairs of samples. Samples with a correlation lower than 0.6 are removed
prep_data <- TCGAanalyze_Preprocessing(object = rse_coding, cor.cut = 0.6, datatype = 'HTSeq - Counts')

# Filter TCGA samples according to tumor purity
purity <- TCGAtumor_purity(colnames(prep_data), 0, 0, 0, 0, 0.6)


# Create matrix containing pure samples as well as samples with no purity info (normal)
puri_data <- prep_data[,colnames(prep_data) %in% unlist(purity)] 


# Normalisation adjusts global properties of measurements for individual samples so that they can be more appropriately compared. In this way, it is ensured that the expression distribution of samples is similar across the entire experiment
norm_data <- TCGAanalyze_Normalization(tabDF =  puri_data,
                                       geneInfo = geneInfoHT,
                                       method = 'gcContent')
norm_data_id <- as.matrix(norm_data)

# Swop ensembl_gene_id to external_gene_name
rownames(norm_data_id) <- rowData(rse_coding)[match(rownames(norm_data_id), rowData(rse_coding)[,'ensembl_gene_id']),'external_gene_name'] 

# Test if candidate genes in matrix
#candidate_test(norm_data_id, candidate_genes)
#candidate_test(norm_data_id, BCL2_genes)


# Filter out genes with very low counts across all libraries, as they provide more noise than evidence, e.g. for differential expression. Genes that are not expressed at a biologically meaningful level in any condition should be discarded to reduce the subset of genes to those that are of interest, and to minimise the number of tests carried out downstream.
filt_data <- TCGAanalyze_Filtering(tabDF = norm_data_id,
                                   method = 'quantile',
                                   qnt.cut =  0.25)

#candidate_test(filt_data, candidate_genes)  # Test if candidate genes in matrix
#candidate_test(norm_data_id, BCL2_genes)

# create a data frame of TCGA barcode data elements
barcode_elements <- get_IDs(filt_data)

# get subtype information
barcode_elements$subtype <- data_subt[match(barcode_elements$patient, data_subt$patient), 'PAM50.mRNA']
barcode_elements$subtype <- as.factor(as.character(ifelse(barcode_elements$condition == 'normal', as.character('Normal'), as.character(barcode_elements$subtype))))

# remove samples where subtype is not known and Normal-like samples
barcode_elements <- barcode_elements[-which(is.na(barcode_elements$subtype)), ] 
barcode_elements <- barcode_elements[-which(barcode_elements$subtype=='Normal-like'), ]

# add information about the year the sample was taken
names(clinic_df)[names(clinic_df) == 'bcr_patient_barcode'] <- 'patient'
clinic_df$age_at_diag_year <- floor(clinic_df$age_at_diagnosis/365)
clinic_df$diag_year <- clinic_df$age_at_diag_year + 
  clinic_df$year_of_birth
diag_year_df <- clinic_df[, c('patient', 'diag_year')]  # add year the sample was taken

barcode_elements$year <- diag_year_df[match(barcode_elements$patient, diag_year_df$patient), 'diag_year']

barcode_elements$year <- as.factor(as.character(ifelse(is.na(barcode_elements$year), as.character('NA'), as.character(barcode_elements$year))))

# remove samples from count matrix that did not contain subtype information
filt_data <- filt_data[ ,match(barcode_elements$barcode, colnames(filt_data))]

condition <- as.factor(as.character(barcode_elements$condition))
subtype <- as.factor(as.character(barcode_elements$subtype))
year <- as.factor(as.character(barcode_elements$year))
tss <- as.factor(as.character(barcode_elements$tss))
plate <- as.factor(as.character(barcode_elements$plate))

normal <- which(barcode_elements$condition=='normal')
cancer <- which(barcode_elements$condition=='cancer')

#boxplot(filt_data['BCL2',normal],filt_data['BCL2', cancer], names = c('normal', 'cancer')) # test expressions

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### EXPLORATORY ANALYSIS ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# log2 transform data before plotting
log_data <- log2(filt_data + 1)
source('src/expl_plots.R')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### DIFFERENTIAL EXPRESSION ANALYSIS CONDITIONS CANCER VERSUS NORMAL ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create design matrix 

design_matrix <- model.matrix(~0 + condition + tss)

colnames(design_matrix)[1:2] <- c('cancer', 'normal')

contr_matrix <- makeContrasts('CancervsNormal' = cancer - normal, levels = colnames(design_matrix))  # Construct contrasts

pdf(file = 'figs/mean-variance-voom.pdf', width = 11)
par(mfrow=c(1,2))
voom_data <- voom(counts = filt_data, design = design_matrix, plot = TRUE)  # Removing heteroscedascity from count data with voom - Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights

lm_fit <- lmFit(object = voom_data, design = design_matrix)  # Fit linear model for each gene

contr_fit <- contrasts.fit(fit = lm_fit, contrasts = contr_matrix)  # Compute estimated coefficients and standard errors for a given set of contrasts

ebayes_fit <- eBayes(contr_fit) # Empirical Bayes Statistics

plotSA(ebayes_fit)  # Check the mean-variance relationship of the expression data, after voom.
title(main = 'Transformed: Mean-variance trend')
dev.off()

# Identify which genes are significantly differentially expressed
tt <- topTable(ebayes_fit, coef=1, adjust.method ='fdr', number=nrow(ebayes_fit))
# Filter out up - and downregulated genes on the basis of a logFC of 1 and a FDR of 0.05
up <- tt[tt$logFC >= 1 & tt$adj.P.Val < 0.05, ]
down <- tt[tt$logFC <= -1 & tt$adj.P.Val < 0.05, ]
up$dir <- rep('up', nrow(up))
down$dir <- rep('down', nrow(down))
de_genes_list <- list(up, down)
de_genes <- rbind(up, down)

cancer_vs_normal_cand <- de_genes[rownames(de_genes) %in% candidate_genes, ] 


df_cancer_vs_normal_cand <- cbind(names = rownames(cancer_vs_normal_cand), cancer_vs_normal_cand)
write.table(df_cancer_vs_normal_cand, file = 'output/DE.cand.cancer.vs.normal.txt', sep = '\t', row.names = FALSE)


#candidate_test(de_genes, candidate_genes)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### GENE CO-EXPRESSION NETWORK ANALYSIS #### 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

source('src/co_expression_analysis.R')


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### DIFFERENTIAL EXPRESSION ANALYSIS PAM50 SUBTYPES AND NORMAL ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create design matrix 
design_matrix <- model.matrix(~0 + subtype + tss)
colnames(design_matrix)[1:5] <- c('Basal', 'HER2', 'LumA', 'LumB', 'Normal')
head(design_matrix)
# Construct contrasts
contr_matrix <- makeContrasts(
  BasalvsHER2 = Basal - HER2,
  BasalvsLumA = Basal - LumA,
  BasalvsLumB = Basal - LumB,
  BasalvsNormal = Basal - Normal,
  HER2vsLumA = HER2 - LumA,
  HER2vsLumB = HER2 - LumB,
  HER2vsNormal = HER2 - Normal,
  LumAvsLumB = LumA - LumB,
  LumAvsNormal = LumA - Normal,
  LumBvsNormal = LumB - Normal,
  levels = colnames(design_matrix))



# Removing heteroscedascity from count data with voom - Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights
voom_data <- voom(counts = filt_data, design = design_matrix, plot = TRUE)  

lm_fit <- lmFit(object = voom_data, design = design_matrix)  # Fit linear model for each gene

contr_fit <- contrasts.fit(fit = lm_fit, contrasts = contr_matrix)  # Compute estimated coefficients and standard errors for a given set of contrasts

ebayes_fit <- eBayes(contr_fit) # Empirical Bayes Statistics

summary(decideTests(ebayes_fit))
plotSA(ebayes_fit)  # Check the mean-variance relationship of the expression data, after voom.
title(main = 'Voom transformed: Mean-variance trend')
dev.off()


# Identify which genes and candidate genes are significantly differentially expressed

# DE genes for Basal-like versus HER2-enriched
basal_vs_her2 <- DE_genes(ebayes_object = ebayes_fit, contrast = 1,
                          adjust_method = 'fdr', LFC = 1, FDR = 0.05)
basal_vs_her2_cand <- basal_vs_her2[rownames(basal_vs_her2) %in% candidate_genes, ]
#dim(basal_vs_her2_cand)
df_basal_vs_her2_cand <- cbind(names = rownames(basal_vs_her2_cand), basal_vs_her2_cand)
write.table(df_basal_vs_her2_cand, file = 'output/DE.cand.basal.vs.her2.txt', sep = '\t', row.names = FALSE)


# DE genes for Basal-like versus Luminal A
basal_vs_lumA <- DE_genes(ebayes_object = ebayes_fit, contrast = 2,
                          adjust_method = 'fdr', LFC = 1, FDR = 0.05)
basal_vs_lumA_cand <- basal_vs_lumA[rownames(basal_vs_lumA) %in% candidate_genes, ]
#dim(basal_vs_lumA_cand)
df_basal_vs_lumA_cand <- cbind(names = rownames(basal_vs_lumA_cand), basal_vs_lumA_cand)
write.table(df_basal_vs_lumA_cand, file = 'output/DE.cand.basal.vs.lumA.txt', sep = '\t', row.names = FALSE)

# DE genes for Basal-like versus Luminal B
basal_vs_lumB <- DE_genes(ebayes_object = ebayes_fit, contrast = 3,
                          adjust_method = 'fdr', LFC = 1, FDR = 0.05)
basal_vs_lumB_cand <- basal_vs_lumB[rownames(basal_vs_lumB) %in% candidate_genes, ]
#dim(basal_vs_lumB_cand)
df_basal_vs_lumB_cand <- cbind(names = rownames(basal_vs_lumB_cand), basal_vs_lumB_cand)
write.table(df_basal_vs_lumB_cand, file = 'output/DE.cand.basal.vs.lumB.txt', sep = '\t', row.names = FALSE)

# DE genes for Basal-like versus Normal
basal_vs_norm <- DE_genes(ebayes_object = ebayes_fit, contrast = 4,
                          adjust_method = 'fdr', LFC = 1, FDR = 0.05)
basal_vs_norm_cand <- basal_vs_norm[rownames(basal_vs_norm) %in% candidate_genes, ]
#dim(basal_vs_norm_cand)
df_basal_vs_norm_cand <- cbind(names = rownames(basal_vs_norm_cand), basal_vs_norm_cand)
write.table(df_basal_vs_norm_cand, file = 'output/DE.cand.basal.vs.norm.txt', sep = '\t', row.names = FALSE)


# DE genes for HER2-enriched versus Luminal A
her2_vs_lumA <- DE_genes(ebayes_object = ebayes_fit, contrast = 5,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
her2_vs_lumA_cand <- her2_vs_lumA[rownames(her2_vs_lumA) %in% candidate_genes, ]
#dim(her2_vs_lumA_cand)
df_her2_vs_lumA_cand <- cbind(names = rownames(her2_vs_lumA_cand), her2_vs_lumA_cand)
write.table(df_her2_vs_lumA_cand, file = 'output/DE.cand.her2.vs.lumA.txt', sep = '\t', row.names = FALSE)

# DE genes for HER2-enriched versus Luminal B
her2_vs_lumB <- DE_genes(ebayes_object = ebayes_fit, contrast = 6,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
her2_vs_lumB_cand <- her2_vs_lumB[rownames(her2_vs_lumB) %in% candidate_genes, ]
#dim(her2_vs_lumB_cand)
df_her2_vs_lumB_cand <- cbind(names = rownames(her2_vs_lumB_cand), her2_vs_lumB_cand)
write.table(df_her2_vs_lumB_cand, file = 'output/DE.cand.her2.vs.lumB.txt', sep = '\t', row.names = FALSE)

# DE genes for HER2-enriched versus Normal
her2_vs_norm <- DE_genes(ebayes_object = ebayes_fit, contrast = 7,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
her2_vs_norm_cand <- her2_vs_norm[rownames(her2_vs_norm) %in% candidate_genes, ]
#dim(her2_vs_norm_cand)
df_her2_vs_norm_cand <- cbind(names = rownames(her2_vs_norm_cand), her2_vs_norm_cand)
write.table(df_her2_vs_norm_cand, file = 'output/DE.cand.her2.vs.norm.txt', sep = '\t', row.names = FALSE)

# DE genes for Luminal A versus Luminal B
lumA_vs_lumB <- DE_genes(ebayes_object = ebayes_fit, contrast = 8,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
lumA_vs_lumB_cand <- lumA_vs_lumB[rownames(lumA_vs_lumB) %in% candidate_genes, ]
#dim(lumA_vs_lumB_cand)
df_lumA_vs_lumB_cand <- cbind(names = rownames(lumA_vs_lumB_cand),lumA_vs_lumB_cand)
write.table(df_lumA_vs_lumB_cand, file = 'output/DE.cand.lumA.vs.lumB.txt', sep = '\t', row.names = FALSE)

# DE genes for Luminal A versus Normal
lumA_vs_norm <- DE_genes(ebayes_object = ebayes_fit, contrast = 9,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
lumA_vs_norm_cand <- lumA_vs_norm[rownames(lumA_vs_norm) %in% candidate_genes, ]
#dim(lumA_vs_norm_cand)
df_lumA_vs_norm_cand <- cbind(names = rownames(lumA_vs_norm_cand),lumA_vs_norm_cand)
write.table(df_lumA_vs_norm_cand, file = 'output/DE.cand.lumA.vs.norm.txt', sep = '\t', row.names = FALSE)

# DE genes for Luminal B versus Normal
lumB_vs_norm <- DE_genes(ebayes_object = ebayes_fit, contrast = 10,
                         adjust_method = 'fdr', LFC = 1, FDR = 0.05)
lumB_vs_norm_cand <- lumB_vs_norm[rownames(lumB_vs_norm) %in% candidate_genes, ]
#dim(lumB_vs_norm_cand)
df_lumB_vs_norm_cand <- cbind(names = rownames(lumB_vs_norm_cand),lumB_vs_norm_cand)
write.table(df_lumB_vs_norm_cand, file = 'output/DE.cand.lumB.vs.norm.txt', sep = '\t', row.names = FALSE)

#candidate_test(basal_vs_her2_cand, candidate_genes)  # Test if candidate genes are excluded
#candidate_test(basal_vs_her2_cand, BCL2_genes)  # Test if candidate genes are excluded


# Create table of regulation directions for all DE candidate genes
deg_all <- data.frame(matrix(ncol = 0, nrow = 295))
rownames(deg_all) <- candidate_genes

deg_all$CancervsNormal <- cancer_vs_normal_cand[match(rownames(deg_all), rownames(cancer_vs_normal_cand)), 'dir']
deg_all$BasalvsHER2 <- basal_vs_her2_cand[match(rownames(deg_all), rownames(basal_vs_her2_cand)), 'dir']
deg_all$BasalvsLumA <- basal_vs_lumA_cand[match(rownames(deg_all), rownames(basal_vs_lumA_cand)), 'dir']
deg_all$BasalvsLumB <- basal_vs_lumB_cand[match(rownames(deg_all), rownames(basal_vs_lumB_cand)), 'dir']
deg_all$BasalvsNormal <- basal_vs_norm_cand[match(rownames(deg_all), rownames(basal_vs_norm_cand)), 'dir']
deg_all$HER2vsLumA <- her2_vs_lumA_cand[match(rownames(deg_all), rownames(her2_vs_lumA_cand)), 'dir']
deg_all$HER2vsLumB <- her2_vs_lumB_cand[match(rownames(deg_all), rownames(her2_vs_lumB_cand)), 'dir']
deg_all$HER2vsNormal <- her2_vs_norm_cand[match(rownames(deg_all), rownames(her2_vs_norm_cand)), 'dir']
deg_all$LumAvsLumB <- lumA_vs_lumB_cand[match(rownames(deg_all), rownames(lumA_vs_lumB_cand)), 'dir']
deg_all$LumAvsNormal <- lumA_vs_norm_cand[match(rownames(deg_all), rownames(lumA_vs_norm_cand)), 'dir']
deg_all$LumBvsNormal <- lumB_vs_norm_cand[match(rownames(deg_all), rownames(lumB_vs_norm_cand)), 'dir']

# # add column of known bh-3 containing id's
# other_bh3 <- read.table(file = 'data/other.known.bh3.txt')
# deg_all$claimedBH3 <- other_bh3[match(rownames(deg_all), other_bh3$V1), 'V1']
# deg_all$claimedBH3 <- ifelse(is.na(deg_all$claimedBH3), NA, 'YES')

#write.table(deg_all, file = 'output/DE.all.txt', sep = '\t')
write.xlsx(deg_all, file = 'output/DE.all.xlsx', rowNames = TRUE)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### PLOTTING RESULTS ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

source('src/results_plot.R')

end_time <- Sys.time()

message(paste0('Analysis completed: ',Sys.Date(), '. Excecution time: ', round(end_time - start_time, 2), ' minutes'))





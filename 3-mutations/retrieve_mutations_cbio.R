#!/usr/bin/env Rscript

#setwd('/Users/simon/master_cbl/20-11-2017_task_2/retrieve_mutations')
setwd(getwd())

# create argument for input file
args = commandArgs(trailingOnly=TRUE)
# test if there is one argument: if not, return an error
if (length(args)==0) {
  stop("A file containing genes of interest must be supplied.\n", call.=FALSE)
}

# Install and load packages if not installed. Othwewise load.
source('https://bioconductor.org/biocLite.R')
list.of.packages <- c('openxlsx', 'cgdsr')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

query_genes <- readLines(args[1])
#query_genes <- c('BCL2A1','HRK', 'GZMB', 'PCNA')
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")

getCancerStudies(mycgds)

# Study: Breast Invasive Carcinoma (TCGA, Provisional)
my_study <- getCancerStudies(mycgds)
my_study <- 'brca_tcga'
case_list <- getCaseLists(mycgds,my_study)[2,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[12,1] # mutations
mutation_data_brca_tcga <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)
my_study <- 'brca_metabric'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[4,1] # mutations
mutation_data_brca_metabric <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast Invasive Carcinoma (British Columbia, Nature 2012)
my_study <- 'brca_bccrc'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[1,1] 
mutation_data_brca_bccrc <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast Invasive Carcinoma (Broad, Nature 2012)
my_study <- 'brca_broad'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[1,1]
mutation_data_brca_broad <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast Invasive Carcinoma (Sanger, Nature 2012)
my_study <- 'brca_sanger'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[1,1]
mutation_data_brca_sanger <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)


# Study: Breast Invasive Carcinoma (TCGA, Cell 2015)
my_study <- 'brca_tcga_pub2015'
case_list <- getCaseLists(mycgds,my_study)[2,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[11,1]
mutation_data_brca_tcga_pup2015 <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast Invasive Carcinoma (TCGA, Nature 2012)
my_study <- 'brca_tcga_pub'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[8,1]
mutation_data_brca_tcga_pup <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Breast cancer patient xenografts (British Columbia, Nature 2014)
my_study <- 'brca_bccrc_xenograft_2014'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[3,1]
mutation_data_brca_bccrc_xenograft_2014 <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: The Metastatic Breast Cancer Project (Provisional, October 2017)
my_study <- 'brca_mbcproject_wagle_2017'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[2,1]
mutation_data_brca_mbcproject_wagle_2017 <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)

# Study: Adenoid Cystic Carcinoma of the Breast (MSKCC, J Pathol. 2015)
my_study <- 'acbc_mskcc_2015'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[2,1]
mutation_data_acbc_mskcc_2015 <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)


# Study: Mutational profiles of metastatic breast cancer (France, 2016)
my_study <- 'brca_igr_2015'
case_list <- getCaseLists(mycgds,my_study)[1,1]
genetic_profile <- getGeneticProfiles(mycgds,my_study)[2,1]
mutation_data_brca_igr_20155 <- getMutationData(mycgds, caseList = case_list, geneticProfile = genetic_profile,  genes = query_genes)


all_mutations <- rbind(mutation_data_acbc_mskcc_2015, mutation_data_brca_bccrc,mutation_data_brca_bccrc_xenograft_2014, mutation_data_brca_broad, mutation_data_brca_igr_20155,
                       mutation_data_brca_mbcproject_wagle_2017, mutation_data_brca_metabric, mutation_data_brca_sanger, mutation_data_brca_tcga, mutation_data_brca_tcga_pup,
                       mutation_data_brca_tcga_pup2015)

write.table(all_mutations, file = 'mutations_cbio.txt', sep = '\t', row.names = FALSE)
write.xlsx(all_mutations, file = 'mutations_cbio.xlsx')




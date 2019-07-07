
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### DOWNLOAD DATA #### 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# harmonized
query_harmo <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       workflow.type = 'HTSeq - Counts',
                       sample.type = c('Primary solid Tumor','Solid Tissue Normal'))

GDCdownload(query_harmo)
data_harmo <- GDCprepare(query_harmo,save = TRUE, save.filename = paste0('data/brca_exp_harmonized_',Sys.Date(),'.rda'))

query_clin <- GDCquery_clinic(project = 'TCGA-BRCA', type = 'clinical')
save(query_clin, file = 'data/BRCA/brca_clin.RData')


dim(data_harmo)
length(which(colData(data_harmo)$shortLetterCode =='TP'))
length(which(colData(data_harmo)$shortLetterCode =='NT'))



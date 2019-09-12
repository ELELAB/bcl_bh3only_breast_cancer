library(SummarizedExperiment)
library(TCGAbiolinks)
library(ggplot2)


cancer_types <- c("BRCA",
                  "GBM",
                  "LUAD",
                  "UCEC",
                  "KIRC",
                  "HNSC",
                  "THCA", 
                  "LUSC",
                  "PRAD",
                  "COAD",
                  "STAD",
                  "BLCA",
                  "LIHC",
                  "KIRP",
                  "ESCA",
                  "READ",
                  "KICH",
                  "CHOL"
)

mutname <- c("p.L99R", "p.M75R", "p.Y120C")
pipe <- "mutect2"
all_muts <- data.frame()
for(cancer in cancer_types){
  query.maf <- GDCquery_Maf(cancer, pipelines = pipe)
  mutpos <- which(query.maf$HGVSp_Short %in% mutname)
  query.mut <- query.maf[mutpos,]
  query.mut <- query.mut[,c(1,37,42,73,74,94)]
  query.mut$cancer_type <- rep(cancer, nrow(query.mut))
  all_muts <- rbind(all_muts, query.mut)
                               
}
write.csv(all_muts, file = "all_cancer_mut.csv")

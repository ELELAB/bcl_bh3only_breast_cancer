library(SummarizedExperiment)
library(TCGAbiolinks)
library(ggplot2)

setwd("8.TCGA_somatic_variations")
source("TCGA_somatic_variations_functions.r")
dir.create("plots", showWarnings = FALSE)
cancer_types <- c("BRCA",
)

geneNames <- c("BCL2A1",
)


all_mut_table <- data.frame()
for(cancer in cancer_types){
  # Find somatic variations using mutect2
  query_maf <- mutect2_tcga_maf(cancer)
  query_gen <- query_maf[which(query_maf$Hugo_Symbol %in% geneNames),]
  
  if(nrow(query_gen) > 0){
    
    # Make a table with impact, polyphen and sift imformation 
    mutations_table <- as.data.frame(rep(cancer, length(geneNames)))
    colnames(mutations_table) <- "Cancer"
    mutations_table$Hugo_Symbol <- as.character(geneNames)
    mutations_table$IMPACT_High <- as.numeric(rep(0, nrow(mutations_table)))
    mutations_table$IMPACT_Moderate <- as.numeric(rep(0, nrow(mutations_table)))
    mutations_table$PolyPhen_PR <- as.numeric(rep(0, nrow(mutations_table)))
    mutations_table$PolyPhen_PO <- as.numeric(rep(0, nrow(mutations_table)))
    mutations_table$SIFT_deliterious <- as.numeric(rep(0, nrow(mutations_table)))
    mutations_table$SIFT_deliterious_low_confidence <- as.numeric(rep(0, nrow(mutations_table)))
    
    for(gene  in query_gen$Hugo_Symbol){
      genpos <- which(query_gen$Hugo_Symbol == gene)
      mutgen_pos <- which(mutations_table$Hugo_Symbol == gene)
      imp_h_len <- as.numeric(length(which(query_gen[genpos,]$IMPACT == "HIGH")))
      if(imp_h_len > 0){
        mutations_table[mutgen_pos,]$IMPACT_High <- as.numeric(imp_h_len)
      } 
      imp_m_len <- length(which(query_gen[genpos,]$IMPACT == "MODERATE"))
      if(imp_m_len > 0){
        mutations_table[mutgen_pos,]$IMPACT_Moderate <- as.numeric(imp_m_len)
      }
      pp_pr_len <- length(grep("probably_damaging", query_gen[genpos,]$PolyPhen))
      if(pp_pr_len > 0){
        mutations_table[mutgen_pos,]$PolyPhen_PR <- as.numeric(pp_pr_len)
      }
      pp_po_len <- length(grep("possibly_damaging", query_gen[genpos,]$PolyPhen))
      if(pp_po_len > 0){
        mutations_table[mutgen_pos,]$PolyPhen_PO <- as.numeric(pp_po_len)
      }
      sif_d_len <- length(grep("deleterious", query_gen[genpos,]$SIFT))
      if(sif_d_len > 0){
        mutations_table[mutgen_pos,]$SIFT_deliterious <- as.numeric(sif_d_len)
      }
      sif_del_len <- length(grep("deleterious_low_confidence", query_gen[genpos,]$SIFT))
      if(sif_del_len > 0){
        mutations_table[mutgen_pos,]$deliterious_low_confidence <- as.numeric(sif_del_len)
      }
      
      
    }
    plot_mut <- melt(mutations_table)
    # Plot results
    #filename_barplot = filename of .pdf barplot file
    filename_barplot <- paste0("plots/",cancer, "_somatic_variations_barplot.pdf")
    ggplot(plot_mut, aes(x=Hugo_Symbol, y=value, fill=variable), cex.axis=1, cex.names= 0.6)+geom_bar(stat='identity', position='dodge')+
      xlab("Genes") + ylab("Number of mutations") + ggtitle(paste0(cancer, " somatic variations"))
    ggsave(filename_barplot)
  }
  all_mut_table <- rbind(all_mut_table, mutations_table)
  write.csv(all_mut_table, file = "tables/impact_polyphen_sift_imformation.csv")

}
pipeline <- c("mutect2", "varscan2", "somaticsniper", "muse")
possposs <- c(1,5,6,7,8,9,35,36,52,53,54,55,56,57,61, 69,73,74,94, 121, 122)
all_muts_pipe <- data.frame()
for(cancer in cancer_types){
  # Find somatic variations using mutect2
  for(pipe in pipeline){
    query_maf <- mutect2_tcga_maf(cancer, pipe)
    query_gen <- query_maf[which(query_maf$Hugo_Symbol %in% geneNames),]
    #query_gen <- query_gen[,possposs]
    query_gen$Cancer <- rep(cancer, nrow(query_gen))
    query_gen$Pipeline <- rep(pipe, nrow(query_gen))
    all_muts_pipe <- rbind(all_muts_pipe, query_gen)

  }
  
}
save(all_muts_pipe, file = "all_mutationspipelines_all_cols.rda")
all_muts_pipe_good <- all_muts_pipe[,possposs]
all_muts_pipe_good <- all_muts_pipe_good[,c(21,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
write.csv(all_muts_pipe_good, file = "all_pipelines_good.csv")

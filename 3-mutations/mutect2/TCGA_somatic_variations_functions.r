library(ggplot2)
library(reshape2)

mutect2_tcga_maf <- function(cancer, pipe){
  dir.create("data", showWarnings = FALSE)
  dir.create("tables", showWarnings = FALSE)
  query.maf <- GDCquery_Maf(cancer, pipelines = pipe, directory = "data/GDCdata")
  #save(query.maf, file = paste0("tables/",cancer, "_somatic_variations.rda"))
  return(query.maf)
}

plot_moderate_high_impact <- function(query.maf, geneNames, colours_genes){
  query.maf <- as.data.frame(query.maf)
  dir.create("plots", showWarnings = FALSE)
  query.maf_genes <- query.maf[which(query.maf$Hugo_Symbol %in% geneNames),]
  if(nrow(query.maf_genes) > 0){
    save(query.maf_genes, file = paste0("tables/",cancer,"_somatic_variations_genes.rda"))
    impacts_to_use <- c("MODERATE", "HIGH")
    gen_imp <- c(which(colnames(query.maf_genes) == "Hugo_Symbol"),which(colnames(query.maf_genes) == "IMPACT"))
    impacts_genes <- query.maf_genes[,gen_imp]
    impact_table <- impacts_genes[which(impacts_genes$IMPACT %in% impacts_to_use),]
    genes <- unique(impact_table$Hugo_Symbol)
    df <- as.data.frame(geneNames)
    colnames(df) <- "genes"
    df$genes <- as.character(df$genes)
    df$number_of_mutations <- c(rep(0, nrow(df)))
    for(gen in df$genes){
      gen_len <- length(which(impact_table$Hugo_Symbol == gen))
      df[which(df$genes == gen),]$number_of_mutations <- gen_len
    }
    genes_true <- geneNames[which(geneNames %in% df$genes)]
    genes_order <- c()
    for(genes in genes_true){
      mut <- which(df$genes == genes)
      genes_order <- c(genes_order, mut)
    }
    df <- df[genes_order,]
    col_order <- which(geneNames %in% df$genes)
    cols <- colours_genes[col_order]
    length_yaxis <- sort(df$number_of_mutations)[nrow(df)]
    if(length_yaxis > 9){
      length_yaxis <- length_yaxis + 5
    }
    if(10 > length_yaxis && length_yaxis > 5){
      length_yaxis <- length_yaxis + 3
    }
    if(length_yaxis <= 5){
      length_yaxis <- length_yaxis + 1
    }
    pdf(paste0("plots/",cancer,"_somatic_variations_barplot.pdf"))
    bp <- barplot(df$number_of_mutations,
                  names.arg = df$genes,
                  ylim = c(0,length_yaxis),
                  cex.axis=1,
                  cex.names=0.7,
                  main = paste(cancer, " somatic variations"),
                  ylab = "Number of somatic variations",
                  col = cols,
                  axes = TRUE)
    dev.off()
  } else {
    cat("\n","No somatic variations was found for your genes in ", cancer, "\n")
  }
  
  
}

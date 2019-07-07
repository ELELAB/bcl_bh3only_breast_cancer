#### FUNCTIONS ####
 
# Test if candidate genes are excluded - print
candidate_test <- function(matrix_like, candidate_vector) {
  for(i in candidate_vector) {
    if (!i %in% rownames(matrix_like)) {
      print(paste0(i,' is excluded!'))
    }
  }
}

# Test if candidate genes are excluded - return vector
get_candidates_excl <- function(matrix_like, candidate_vector) {
  cand <- character()
  for(i in seq_along(candidate_vector)) {
    if (!candidate_vector[i] %in% rownames(matrix_like)) {
      cand <- c(cand, candidate_vector[i])
    }
  }
  return(cand)
} 

# Create adjacency data frame from matrix
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flatten_corr_matrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Function to create MDS plot of condition 
mds_plot_condition <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face='bold'), axis.title=element_text(size=16,face='bold')) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = 'top') + theme(axis.text=element_text(size=16, face='bold')) + theme(axis.text = element_text(colour = 'black'))
}

# Function to create MDS plot of plate (Explore batch effect)
mds_plot_plate <- function(my.data, my.group, my.labels) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face='bold'), axis.title=element_text(size=16,face='bold')) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = 'top') + theme(axis.text=element_text(size=16, face='bold')) + theme(axis.text = element_text(colour = 'black'))
}

# Function to create MDS plot of tss (Explore batch effect)
mds_plot_tss <- function(my.data, my.group, my.labels) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face='bold'), axis.title=element_text(size=16,face='bold')) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = 'top') + theme(axis.text=element_text(size=16, face='bold')) + theme(axis.text = element_text(colour = 'black'))
}
# Function to create MDS plot w. normal confidence ellipses
mds_plot_ellipses <- myMDSplot <- function(my.data, my.group) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  ggplot(data=res, aes(x=M1,y=M2,color=my.group)) + stat_ellipse(type = 'norm') + stat_ellipse(type = 'euclid', geom = 'polygon', level = 2, aes(fill = my.group)) + coord_fixed() + coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank())
}

# This function splits the RangedSummarizedExperiment object and returns a data frame of TCGA data elements (identifiers)
get_IDs <- function(data) {
  # splits column names by '-'
  IDs <- strsplit(c(colnames(data)), '-')
  # for each element in IDs list, combine by row into data frame
  IDs <- ldply(IDs, rbind)
  # create column names
  colnames(IDs) <- c('project', 'tss','participant', 'sample', 'portion', 'plate', 'center')
  cols <- c('project', 'tss', 'participant')
  # add patiet tag
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = '-' )
  barcode <- colnames(data)
  # add barcode to data frame
  IDs <- cbind(IDs, barcode)
  condition <- gsub('11+[[:alpha:]]', 'normal', as.character(IDs$sample))
  condition  <- gsub('01+[[:alpha:]]', 'cancer', condition)
  # add the column condition
  IDs$condition <- condition
  # add rownumber 
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
}
# This function identify which genes are significantly differentially expressed and
# filter out up - and downregulated genes on the basis of a logFC and a FDR
DE_genes <- function(ebayes_object, contrast, adjust_method, LFC, FDR) {
  top <- topTable(ebayes_object, coef=contrast, adjust.method = adjust_method, number=nrow(ebayes_object))
  up <- top[top$logFC >= LFC & top$adj.P.Val < FDR, ]
  down <- top[top$logFC <= -LFC & top$adj.P.Val < FDR, ]
  up$dir <- rep('up', nrow(up))
  down$dir <- rep('down', nrow(down))
  de_list <- list(up, down)
  de_genes <- rbind(up, down)
  return(de_genes)
}

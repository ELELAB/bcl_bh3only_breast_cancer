setwd(getwd())


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### EXPLORATORY ANALYSIS ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# PCA showing the variance along PC1 and PC2 for PAM50 molecular subtypes of breast cancer and normal samples.
pdf(file = 'figs/pca.pam50.subtype.pdf')
pca <- prcomp(t(log_data), scale = TRUE) 
g <- ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, 
              groups = subtype, ellipse = TRUE, circle = TRUE, var.axes = FALSE) + ggtitle('PCA plot of PAM50 cancer subtypes and normal samples')
print(g)
dev.off()

# PCA showing the variance along PC1 and PC2 for cancer and normal samples.
pdf(file = 'figs/pca.cancervsnormal.pdf')
g <- ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, 
         groups = condition, ellipse = TRUE, circle = TRUE, var.axes = FALSE) + ggtitle('PCA plot of cancer and normal samples') + theme(plot.title = element_text(hjust = 0.5))
print(g)
dev.off()

# To assess a possible batch effect we carry out PCA of year, tss and plate
pdf(file = 'figs/pca.year.pdf')
g <- ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, 
         groups = year, ellipse = TRUE, circle = TRUE, var.axes = FALSE) + ggtitle('PCA plot of year') + theme(plot.title = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(file = 'figs/pca.tss.pdf')
g <- ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, 
         groups = tss, ellipse = TRUE, circle = TRUE, var.axes = FALSE) + ggtitle('PCA plot of tissue source site') + theme(plot.title = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(file = 'figs/pca.plate.pdf')
g <- ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, 
         groups = plate, ellipse = TRUE, circle = TRUE, var.axes = FALSE) + ggtitle('PCA plot of plate') + theme(plot.title = element_text(hjust = 0.5))
print(g)
dev.off()

## Because there are so many samples in our study it becomes very diffucult to spot a possuble clustering trend. Instead of two-dimensional PCA we can explore the variation using One-dimensional PCA, thus reducing the noise of the many samples
# s <- svd(log_data)
# # PC1-4 year
# pdf(file = 'figs/1D.pca.year.pdf')
# tit <- c('PC1-year', 'PC2-year','PC3-year','PC4-year')
# mypar(2,2)
# for(i in 1:4){
#   boxplot(split(s$v[,i],year),las=2,range=0, main = tit[i])
#   stripchart(split(s$v[,i],year),add=TRUE,vertical=TRUE,pch=1,cex=0.5,col=1)
# }
# dev.off()
# # PC1-4 tss
# pdf(file = 'figs/1D.pca.tss.pdf')
# tit <- c('PC1-tss', 'PC2-tss','PC3-tss','PC4-tss')
# mypar(2,2)
# for(i in 1:4){
#   boxplot(split(s$v[,i],tss),las=2,range=0, main = tit[i])
#   stripchart(split(s$v[,i], tss),add=TRUE,vertical=TRUE,pch=1,cex=.5,col=1)
# }
# dev.off()
# # PC1-4 plate
# pdf(file = 'figs/1D.pca.plate.pdf')
# tit <- c('PC1-plate', 'PC2-plate','PC3-plate','PC4-plate')
# mypar(2,2)
# for(i in 1:4){
#   boxplot(split(s$v[,i], plate),las=2,range=0, main = tit[i])
#   stripchart(split(s$v[,i], plate),add=TRUE,vertical=TRUE,pch=1,cex=.5,col=1)
# }
# dev.off()


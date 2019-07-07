setwd(getwd())
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### PLOTTING RESULTS ####
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create volcano plot where significantly differentially expressedcandidate genes are highlighted ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cancer versus normal
TCGAVisualize_volcano(x = tt$logFC,
                      y = tt$adj.P.Val,
                      highlight = rownames(tt[match(rownames(cancer_vs_normal_cand),rownames(tt)), ]),
                      filename = 'figs/cancer.vs.normal.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(tt),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Cancer vs. Normal)',
                      width = 10)
# Basal-like versus HER2-enriched
TCGAVisualize_volcano(x = basal_vs_her2$logFC,
                      y = basal_vs_her2$adj.P.Val,
                      highlight = rownames(basal_vs_her2[match(rownames(basal_vs_her2_cand),rownames(basal_vs_her2)), ]),
                      filename = 'figs/basal.vs.her2.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(basal_vs_her2),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Basal-like vs. HER2-enriched)',
                      width = 10)
# Basal-like versus Luminal A
TCGAVisualize_volcano(x = basal_vs_lumA$logFC,
                      y = basal_vs_lumA$adj.P.Val,
                      highlight = rownames(basal_vs_lumA[match(rownames(basal_vs_lumA_cand),rownames(basal_vs_lumA)), ]),
                      filename = 'figs/basal.vs.lumA.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(basal_vs_lumA),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Basal-like vs. Luminal A)',
                      width = 10)
# Basal-like versus Luminal B
TCGAVisualize_volcano(x = basal_vs_lumB$logFC,
                      y = basal_vs_lumB$adj.P.Val,
                      highlight = rownames(basal_vs_lumB[match(rownames(basal_vs_lumB_cand),rownames(basal_vs_lumB)), ]),
                      filename = 'figs/basal.vs.lumB.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(basal_vs_lumB),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Basal-like vs. Luminal B)',
                      width = 10)
# Basal-like versus Normal
TCGAVisualize_volcano(x = basal_vs_norm$logFC,
                      y = basal_vs_norm$adj.P.Val,
                      highlight = rownames(basal_vs_norm[match(rownames(basal_vs_norm_cand),rownames(basal_vs_norm)), ]),
                      filename = 'figs/basal.vs.normal.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(basal_vs_norm),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Basal-like vs. Normal)',
                      width = 10)
# HER2-enriched versus Luminal A
TCGAVisualize_volcano(x = her2_vs_lumA$logFC,
                      y = her2_vs_lumA$adj.P.Val,
                      highlight = rownames(her2_vs_lumA[match(rownames(her2_vs_lumA_cand),rownames(her2_vs_lumA)), ]),
                      filename = 'figs/her2.vs.lumA.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(her2_vs_lumA),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (HER2-enriched vs. Luminal A)',
                      width = 10)
# HER2-enriched versus Luminal B
TCGAVisualize_volcano(x = her2_vs_lumB$logFC,
                      y = her2_vs_lumB$adj.P.Val,
                      highlight = rownames(her2_vs_lumB[match(rownames(her2_vs_lumB_cand),rownames(her2_vs_lumB)), ]),
                      filename = 'figs/her2.vs.lumB.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(her2_vs_lumB),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (HER2-enriched vs. Luminal B)',
                      width = 10)
# HER2-enriched versus Normal
TCGAVisualize_volcano(x = her2_vs_norm$logFC,
                      y = her2_vs_norm$adj.P.Val,
                      highlight = rownames(her2_vs_norm[match(rownames(her2_vs_norm_cand),rownames(her2_vs_norm)), ]),
                      filename = 'figs/her2.vs.normal.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(her2_vs_norm),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (HER2-enriched vs. Normal)',
                      width = 10)
# Luminal A versus Luminal B
TCGAVisualize_volcano(x = lumA_vs_lumB$logFC,
                      y = lumA_vs_lumB$adj.P.Val,
                      highlight = rownames(lumA_vs_lumB[match(rownames(lumA_vs_lumB_cand),rownames(lumA_vs_lumB)), ]),
                      filename = 'figs/lumA.vs.lumB.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(lumA_vs_lumB),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Luminal A vs. Luminal B)',
                      width = 10)

# Luminal A versus Normal
TCGAVisualize_volcano(x = lumA_vs_norm$logFC,
                      y = lumA_vs_norm$adj.P.Val,
                      highlight = rownames(lumA_vs_norm[match(rownames(lumA_vs_norm_cand),rownames(lumA_vs_norm)), ]),
                      filename = 'figs/lumA.vs.norm.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(lumA_vs_norm),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Luminal A vs. Normal)',
                      width = 10)
# Luminal B versus Normal
TCGAVisualize_volcano(x = lumB_vs_norm$logFC,
                      y = lumB_vs_norm$adj.P.Val,
                      highlight = rownames(lumB_vs_norm[match(rownames(lumB_vs_norm_cand),rownames(lumB_vs_norm)), ]),
                      filename = 'figs/lumB.vs.norm.volcanoexp.pdf',
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(lumB_vs_norm),
                      show.names = 'highlighted',
                      color = c('black','red','darkgreen'),
                      names.size = 2,
                      xlab = ' Gene expression fold change (Log2)',
                      legend = 'State',
                      title = 'Volcano plot (Luminal B vs. Normal)',
                      width = 10)



### Create heatmap of differentially expressed candidate genes between cancer vs. normal comparison ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(file = 'figs/heatmap.cancervsnormal.pdf')
i <- which(rownames(voom_data$E) %in% rownames(cancer_vs_normal_cand))
color = colorRampPalette(rev(brewer.pal(n = 9, name ='RdYlBu')))
col_condition <- ifelse(barcode_elements$condition=='normal','blue','red')
heatmap.2(voom_data$E[i, ],col = color, labCol = condition, trace = 'none',
          density.info = 'none', dendrogram = 'column', hclustfun = function(d) hclust(d, method = 'ward.D2'),
          ColSideColors = col_condition, scale = 'row')
legend('left' ,  
       legend = c('NORMAL', 'CANCER'),
       col = c('blue', 'red'), 
       border = FALSE, bty = 'n', cex = 0.7, lty = 1, lwd = 3)
dev.off()

### Create heatmap of differentially expressed candidate genes between subtype comparisons ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

deg_cand <- deg_all[rowSums(is.na(deg_all)) != ncol(deg_all),]

i <- which(rownames(voom_data$E) %in% rownames(deg_cand))

cols <- brewer.pal(nlevels(subtype), 'Set3')
cols <- factor(subtype, labels = cols)

pdf('figs/heatmap.subtypes.pdf', paper='a4', width=8, height=8)
heatmap.2(voom_data$E[i, ], col = rev(heat.colors(50)), trace = 'none',
          density.info = 'none', dendrogram = 'column', hclustfun = function(d) hclust(d, method = 'ward.D2'),
          ColSideColors = as.character(cols), scale = 'row', labCol = '', main = 'Candidate gene cluestering', cexRow = 0.8)
legend('left' ,  
       legend = c('Basal-like', 'HER2', 'Luminal A', 'Luminal B', 'Normal'),
       col = c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3'), 
       border = FALSE, bty = 'n', cex = 0.7, lty = 1, lwd = 3)
dev.off()


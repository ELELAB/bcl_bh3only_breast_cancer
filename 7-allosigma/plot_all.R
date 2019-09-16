############################################
library(ggplot2)
library(viridis)
library(reshape2)
library(tidyverse)
library(plyr)
################################################################################################################
################################################################################################################
allosigma = read.table("5UUL_model3.txt", header = FALSE, fill = TRUE)
allosigmax <- tibble::rowid_to_column(allosigma, "res.num")
allosigmaxx <- allosigmax[1:151,c(1,3,5)]
allosigmaxxx <- rename(allosigmaxx, c("res.num"="res.num", "V2"="5UUL", "V4"="CABSflex conformation"))
allosigmaxxxx <- melt(allosigmaxxx, id.vars="res.num")
d <- rename(allosigmaxxxx, c("res.num"="res.num", "variable"="structure", "value"="Allo_ener"))
setEPS()
postscript("AlloSIGMA_5UUL_CABSflex.model.eps")
ggplot(d, aes(x=res.num, y=Allo_ener, group=structure)) +
  geom_line(aes(color=structure)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(y="Allosteric free energy (kcal/mol)") +
  theme_classic()
dev.off()
p <- ggplot(d, aes(x=res.num, y=Allo_ener, group=structure)) +
  geom_line(aes(color=structure)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(y="Allosteric free energy (kcal/mol)") +
  theme_classic() 
p2 <- p + geom_hline(yintercept=0, linetype="dashed") + scale_x_continuous(limits = c(1, 151), expand = c(0,0), breaks=c(1, 25, 50, 75, 100, 125, 151))
ggsave("test.2.pdf",p2, width=10, height=3, units="in", scale=1)



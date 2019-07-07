# Analyses of the Cancer Genome Atlas (TCGA) RNASeq BRCA dataset

For this study we aggregated, pre-processed and normalized RNA-Seq breast (BRCA) cancer data from The Cancer Genome Atlas (TCGA), using a new version of the R package [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html). The data in this study is accessible through the [NCI Genomic Data Commons](https://gdc.cancer.gov). The GDC Data Portal provides access to the subset of TCGA data that has been harmonized against GRCh38 (hg38) using GDC Bioinformatics Pipeline which provides methods to the standardization of biospecimen and clinical data, the re-alignment of DNA and RNA sequence data against a common reference genome build GRCh38, and the generation of derived data. 
We carried out differential expression analysis based on a solid statistical framework implemented in the R package [Limma-voom](http://bioconductor.org/packages/release/bioc/html/limma.html) to reveal which Bcl-2 family members and their interaction partners are dysregulated in a given molecular cancer subtype. Subsequently we carried out a gene co-expression network analysis of the candidate genes. We also retrieve known mutations of the anti-apoptotoc BCL2 genes from the [Cbioportal](http://www.cbioportal.org).

The [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) dataset used in this analysis was retreived 2017-11-24. 

**To accomplish this we developed the R script ``analysis.R``. The script is a complete pipeline encompassing the following steps:**

1. Retrieve the data
2. Pre-process, normalize and filter data
3. Explore the data visually - allow exploration of unwanted noise and hidden artifacts such as batch effects.
4. Differential expression analysis
5. Export and plot differently expressed candidate genes
6. Gene co-expression network analysis of the candidate genes

## To reproduce this analysis

### Prerequisites

```
R version 3.3.1 or higher
Rstudio version 1.1.383 or higher        
Bioconductor version 3.6 or higher
```
Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages for manual installation:

CRAN:
```
ggplot2
RColorBrewer
gplots
igraph
openxlsx
heatmap.2
devtools
psych
dplyr
cgdsr
```
Bioconductor:
```
SummarizedExperiment
limma
TCGAbiolinks
biomaRt
```
GitHub:
```
ggbiplot
```

### Running the analysis

```
$ cd bcl2_bh3_breast_cancer/tcga_brca_dea/
$ chmod +x analysis.R
$ nohup ./analysis.R &
```

**Running analysis.R will produce the following files:**
```
tcga_brca_dea/
├── analysis.R
├── data
│   ├── brcaExp_harmonized_2017-11-24.rda
│   ├── brca_clin.RData
│   └── candidate.genes.txt
├── figs
│   ├── basal.vs.her2.volcanoexp.pdf - Volcano plot where significantly differentially expressed candidate genes are highlighted
│   ├── basal.vs.lumA.volcanoexp.pdf - *
│   ├── basal.vs.lumB.volcanoexp.pdf - *
│   ├── basal.vs.normal.volcanoexp.pdf - *
│   ├── cancer.vs.normal.volcanoexp.pdf - *
│   ├── co.expr.dendro.pdf - Dendrogram of clustered communities from co-expression network
│   ├── co.expr.net.pdf - Co-expression network of candidate genes (cutoff: 0.6)
│   ├── co.expr.subnet.pdf - Co-expression subnetwork of candidate genes
│   ├── heatmap.subtypes.pdf - Heatmap of differentially expressed candidate genes - subtype comparison
│   ├── heatmap.cancervsnormal.pdf - Heatmap of differentially expressed candidate genes - cancer versus normal comparison
│   ├── her2.vs.lumA.volcanoexp.pdf - Volcano plot where significantly differentially expressed candidate genes are highlighted
│   ├── her2.vs.lumB.volcanoexp.pdf - *
│   ├── her2.vs.normal.volcanoexp.pdf - *
│   ├── lumA.vs.lumB.volcanoexp.pdf - *
│   ├── lumA.vs.norm.volcanoexp.pdf - *
│   ├── lumB.vs.norm.volcanoexp.pdf - *
│   ├── mean-variance-voom.pdf - Mean - variance relationship before and after voom transformation
│   ├── pca.cancervsnormal.pdf - PCA plot showing the variance along PC1 and PC2 for transcriptional profiles of cancer and normal samples
│   ├── pca.pam50.subtype.pdf - PCA plot showing the variance along PC1 and PC2 for transcriptional profiles of PAM50 molecular subtypes and normal samples
│   ├── pca.plate.pdf - PCA plot of transciptional profiles. Allows for exploration of technical variation deriving from different sequencing plates
│   ├── pca.tss.pdf - PCA plot of transciptional profiles. Allows for exploration of technical variation deriving from different tissue source sites (centers)
│   └── pca.year.pdf PCA plot of transciptional profiles. Allows for exploration of technical variation deriving from different year sample was taken
├── output
│   ├── DE.all.xlsx - Direction of regulation of differentially expressed candidate genes of all comparisons
│   ├── DE.cand.basal.vs.her2.txt - Differentially expressed candidate genes between each comparison (including 'adj.P.Val' and 'logFC')
│   ├── DE.cand.basal.vs.lumA.txt - *
│   ├── DE.cand.basal.vs.lumB.txt - *
│   ├── DE.cand.basal.vs.norm.txt - *
│   ├── DE.cand.cancer.vs.normal.txt - *
│   ├── DE.cand.her2.vs.lumA.txt - *
│   ├── DE.cand.her2.vs.lumB.txt - *
│   ├── DE.cand.her2.vs.norm.txt - *
│   ├── DE.cand.lumA.vs.lumB.txt - *
│   ├── DE.cand.lumA.vs.norm.txt - *
│   ├── DE.cand.lumB.vs.norm.txt - *
│   ├── mutations.antiapop.txt - Known mutations retrieved from cbioportal
│   
└── src
    ├── co_expression_analysis.R
    ├── expl_plots.R
    ├── functions.R
    ├── results_plot.R
  

4 directories, 45 files
```

## Session info

```
R version 3.4.3 (2017-11-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.5 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.0
LAPACK: /usr/lib/lapack/liblapack.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets 
 [9] methods   base     

other attached packages:
 [1] ggbiplot_0.55              scales_0.5.0               plyr_1.8.4                
 [4] cgdsr_1.2.10               dplyr_0.7.4                psych_1.7.8               
 [7] biomaRt_2.34.2             devtools_1.13.4            openxlsx_4.0.17           
[10] igraph_1.1.2               TCGAbiolinks_2.7.21        limma_3.34.6              
[13] SummarizedExperiment_1.8.1 DelayedArray_0.4.1         matrixStats_0.53.0        
[16] Biobase_2.38.0             GenomicRanges_1.30.1       GenomeInfoDb_1.14.0       
[19] IRanges_2.12.0             S4Vectors_0.16.0           BiocGenerics_0.24.0       
[22] gplots_3.0.1               RColorBrewer_1.1-2         ggplot2_2.2.1             
```






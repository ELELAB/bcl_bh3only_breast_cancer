### TCGA ##
The script in this folder will allow to retrieve the missense mutations of interest for this study across different tcga cancer studies. In the output files only the entries for BCL2A1 gene should be considered.

The data have been retrieved by using the TCGAbiolinks function GDCquery_Maf with the mutect2pipeline. The pipeline compares tumour to normal samples and a ‘pool of normal’ to find somatic variations (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/  and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5302158/ ).

## For the Gene Set Enrichment Analysis (GSEA)

This analysis is divided into the Functional Comparison and the Taxonomic Comparison
The Required Input Data
* PICRUST Data (.txt): DESeq Comparison SPT and BPT for imputed metagenome; gene names in rows and log fold change and p value in columns.
* WGS Functional Data (.txt): DESeq Comparison SPT and BPT for Whole Genome Sequencing gene names in rows and log fold change and p value in columns.
* RNA Functional Data (.txt): DESeq Comparison SPT and BPT for RNA Metatranscriptome Sequencing gene names in rows and log fold change and p value in columns.
* 16S Data (.txt): DESeq Comparison SPT and BPT for 16S rRNA gene sequencing; Taxa names in rows and log fold change and p value in columns.
* WGS Taxonomic Data (.txt): DESeq Comparison SPT and BPT for Whole Genome Sequencing Taxa names in rows and log fold change and p value in columns.
* RNA Taxonomic Data (.txt): DESeq Comparison SPT and BPT for RNA Metatranscriptome Sequencing Taxa names in rows and log fold change and p value in columns.

1. Open R
1. Follow Code Combined_Omics_Figure_GIT.r to create The Following Figures
    * #### Figure 2B, Figure 2D, Figure 2F, Figure 3A
1. At the start of the code you will need to load packages required for analysis, for example:
```
library(eulerr)
```
* If you don't have this package you will need to install it, for example:
```
install.packages("eulerr")
```

## For the MultiOmic Comparison
The Required Input Data
* PICRUST comparison (.txt): DESeq Comparison SPT and BPT for 16S rRNA gene sequencing gene names in rows and log fold change and p value in columns.
* WGS comparison (.txt): DESeq Comparison SPT and BPT for Whole Genome Sequencing gene names in rows and log fold change and p value in columns.
* RNA comparison (.txt): DESeq Comparison SPT and BPT for RNA Metatranscriptome Sequencing gene names in rows and log fold change and p value in columns.

1. Open R
1. Follow Code Combined_Omics_Figure_GIT.r to create The Following Figure
  * #### Figure 3B
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(DESeq2)
```
* If you don't have this package you will need to install it, for example:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

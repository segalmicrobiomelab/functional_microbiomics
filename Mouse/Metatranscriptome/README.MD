## For the Mouse RNA Metatranscriptome Sequencing
The Required Input Data
* Gene Abundance (.RData):An abundance table with gene names in rows and samples in columns.
* MetaData (.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code MultiOmics_Metatranscriptome_Mouse_GIT.r to create The Following Figures
  * #### Figure 6A, SFigure 6E, Figure 6F, Figure 6G
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

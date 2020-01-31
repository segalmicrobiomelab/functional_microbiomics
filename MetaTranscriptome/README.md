## For the Human RNA Metatranscriptome Sequencing
The Required Input Data
* Gene Abundance (.txt):An abundance table with KO names in rows and samples in columns.
* MetaData (.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code MultiOmics_Metatranscriptome_GIT.r to create The Following Figures
    * #### Figure 1H, SFigure 1B, Figure 2E, Figure 5D, Figure 5E, SFigure 4B, SFigure 7B
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

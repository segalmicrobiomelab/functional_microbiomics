# This analysis is first done for Functional Comparisons and then Taxonomic Comparisons

## For the Functional Comparisons
The Required Input Data
* PICRUST (.RData):Contains the count tabe from the PICRUST Analysis.
* Metatranscriptome (.RData): Contains the count table for only the KO annotated counts for the Metatranscriptome
* Metagenome (.RData): Contains the count table for only the KO annotated counts for the Metagenome

1. Open R
1. Follow Code Multiomics_PROCRUSTES_GIT.r to create the following Figures:
    * #### SFigure 3, SFigure 4
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

## For the Taxonomic Comparisons
The Required Input Data
* 16S (.RData):Contains the count tabe from the 16S rRNA Analysis.
* Metatranscriptome (.RData): Contains the count table for only the TAXA annotated counts for the Metatranscriptome
* Metagenome (.RData): Contains the count table for only the TAXA annotated counts for the Metagenome

1. Open R
1. Follow Code 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")

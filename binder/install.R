# Just need ggplot2 for example, but allows for further play
# install.packages('tidyverse')
install.packages("dplyr")
install.packages("tidyr")
install.packages("BiocManager")
install.packages("vegan")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("pheatmap")
install.packages("ggrepel")
install.packages("scales")
install.packages("data.table")
install.packages("fBasics")
install.packages("forcats")
install.packages("maptools")

#Install Bioconductor Packages
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("limma")
BiocManager::install("KEGGREST")
BiocManager::install("pathview")
BiocManager::install("phyloseq")

#Install Packages That Required BioConductor
install.packages("omu")
install.packages("pathfindR")

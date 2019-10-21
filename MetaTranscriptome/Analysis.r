#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library(SpiecEasi)
library(vegan)

#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
        axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.01

#Load Meta Data
coldata <- read.delim2("Map.Metatrancriptome.COPD.SmNV.b3.txt", sep="\t")

#Remove Sample with 0 Reads for all Genes
coldata <- coldata[coldata$SampleID!="SmNV.0039.BAL.L.untouched",]

#Order Meta Data by SampleId
coldata <- coldata[order(coldata$SampleID),]

#load Count Data
mycounts <-read.delim2("KEGG_gene_table.txt", sep="\t", row.names=1)

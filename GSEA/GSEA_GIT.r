#Load Libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(eulerr)

#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#//////////////////////////////////////////////////////////////////
#//////////////////////FUNCTIONAL DATA/////////////////////////////
#//////////////////////////////////////////////////////////////////


#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>Metatrans DATA
#//////////////////////////////////////////////////////////////////

#Load Data
x <- read.table("Metatrans_BAL_KEGG_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)

#Create a aranks table based on KO_Number
head(x)
res1 <- x %>% 
   dplyr::select(KO_Number, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(KO_Number) %>% 
   summarize(stat=mean(as.numeric(stat)))
res1

#Deframe the ranks
ranks <- deframe(res1)
head(ranks, 20)


#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>Metagenome DATA
#//////////////////////////////////////////////////////////////////
#Load Data
a<-read.table("Metagenome_BAL_KEGG_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)
#remove any NAs
a <- a[!is.na(a$KO_Subclass_2),]

res2 <- a %>%
         group_by(KO_Subclass_2) %>% 
         summarise(KO_Number = list(unique(KO_Number)))

#Creat a List of Pathways with associated KEGG IDs
pathways.hallmark<- lapply(split(x = res2$KO_Number, f = res2$KO_Subclass_2), unlist)

#Check that the list looks good
pathways.hallmark %>% 
   head() %>% 
   lapply(head)

#Make sure it is a list
class(pathways.hallmark)

#Run GSEA Analysis Comparing Metagenaome and Metatranscriptome
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #1051
#Count Number of KOs in Metatranscriptome Data
length(unique(x[["KO_Number"]])) #6450
#Count Number of KOs in Metagenome Data
length(unique(a[["KO_Number"]])) #1740

#Metatranscriptome (6450) is Number minus the overalp and Metagenome (115) is Number minus the overlap
fit <- euler(c(Metatranscriptome = 5399, Metagenome = 689, "Metatranscriptome&Metagenome" = 1051)) 

#---------------
#----Figure 3A
#---------------
pdf("Metatranscriptome_vs_Metagenome_KEGG_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#7A81FF","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()



#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>PICRUST DATA
#//////////////////////////////////////////////////////////////////

#Load Data
s <- read.table("KEGG_PICRUST_BAL.txt", header=T, sep="\t", row.names=1)

#Create Rank of KOs
head(s)
res3 <- s %>% 
   dplyr::select(KO_Number, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(KO_Number) %>% 
   summarize(stat=mean(as.numeric(stat)))
res3

#DEframe Ranks
ranks3 <- deframe(res3)
head(ranks3, 20)

#Run GSEA Analysis Comparing PICRUST to Metagenome
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks3, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #1025
#Count Number of KOs in 16S Data
length(unique(s[["KO_Number"]])) #6296
#Count Number of KOs in Metagenome Data
length(unique(a[["KO_Number"]])) #1740

#Metatranscriptome (6296) is Number minus the overalp and Metagenome (115) is Number minus the overlap
fit <- euler(c("16S" = 5271, Metagenome = 715, "16S&Metagenome" = 1025)) 

#---------------
#----Figure 3A
#---------------
pdf("16S_vs_Metagenome_KEGG_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#FF2F92","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()



#remove any NAs
s <- s[!is.na(s$KO_Subclass_2),]

res4 <- s %>%
         group_by(KO_Subclass_2) %>% 
         summarise(KO_Number = list(unique(KO_Number)))

#Creat a List of Pathways with associated KEGG IDs
pathways.hallmark<- lapply(split(x = res4$KO_Number, f = res4$KO_Subclass_2), unlist)

#Check that the list looks good
pathways.hallmark %>% 
   head() %>% 
   lapply(head)

#Make sure it is a list
class(pathways.hallmark)

#Run GSEA Analysis comparing Metatranscriptome and PICRUST
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#KeepKOs
keepKOs = gsea$leadingEdge

#Remove Duplicate KOs
keepKOs = keepKOs[!duplicated(keepKOs)]

#If you want all shared KOs from the group with the smallest number of KOs
keepKOs = s$KO_Number
keepKOs = keepKOs[!duplicated(keepKOs)]

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #1053
#Count Number of KOs in 16S Data
length(unique(s[["KO_Number"]])) #6296
#Count Number of KOs in Metagenome Data
length(unique(x[["KO_Number"]])) #6450

#Metatranscriptome (6450) is Number minus the overalp and 16S (6296) is Number minus the overlap
fit <- euler(c("16S" = 5243, Metatranscriptome = 4397, "16S&Metatranscriptome" = 1053)) 

#---------------
#----Figure 3A
#---------------
pdf("16S_vs_Metatranscriptome_KEGG_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#FF2F92","#7A81FF"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()




#//////////////////////////////////////////////////////////////////
#//////////////////////TAXA DATA///////////////////////////////////
#//////////////////////////////////////////////////////////////////

#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>Metatrans DATA
#//////////////////////////////////////////////////////////////////

#Load Data
x <- read.table("Metatranscriptome_BAL_TAXA_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)

#Get Genus Table
library(dplyr)
x = x %>% mutate(Genus=gsub("\\..*","",gs))

#Create Ranks
head(x)
res1 <- x %>% 
   dplyr::select(Genus, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(Genus) %>% 
   summarize(stat=mean(as.numeric(stat)))
res1

#Deframe Ranks
ranks <- deframe(res1)
head(ranks, 20)



#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>Metagenome DATA
#//////////////////////////////////////////////////////////////////
#Load Data
a<-read.table("Metagenome_BAL_TAXA_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)

#Get Genus Table
library(dplyr)
a = a %>% mutate(Genus=gsub("\\..*","",gs))

#make sure numeric variables are numeric
a$log2FoldChange <- as.numeric(a$log2FoldChange)
a$padj <- as.numeric(a$padj)
a <- a[!is.na(a$padj),]
#a <- a[a$padj <= 0.1,] # set pdj threshold
a$fcSign <- sign(a$log2FoldChange)

#Goup by Pathway
res2 <- a %>%
         group_by(KO_Subclass_2) %>% 
         summarise(KO_Number = list(unique(KO_Number)))

#Creat a List of Pathways with associated KEGG IDs
pathways.hallmark<- lapply(split(x = res2$KO_Number, f = res2$KO_Subclass_2), unlist)

#Check that the list looks good
pathways.hallmark %>% 
   head() %>% 
   lapply(head)

#Make sure it is a list
class(pathways.hallmark)

#Run GSEA Analysis comparing Metatranscriptome to Metagenome
fgseaRes <- fgsea(pathways=gmt.file, stats=ranks, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #6
#Count Number of KOs in Metatranscriptome Data
length(unique(x[["Genus"]])) #58
#Count Number of KOs in Metagenome Data
length(unique(a[["Genus"]])) #16

#Metatranscriptome (58) is Number minus the overalp and Metagenome (16) is Number minus the overlap
fit <- euler(c(Metatranscriptome = 52, Metagenome = 10, "Metatranscriptome&Metagenome" = 6)) 

#---------------
#----Figure 2F
#---------------
pdf("Metatranscriptome_vs_Metagenome_TAXA_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#7A81FF","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()


#//////////////////////////////////////////////////////////////////
#------------------------------------------------------>16S DATA
#//////////////////////////////////////////////////////////////////

#Load Data
s <- read.table("16S_BAL_TAXA_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)

#Rank the DAta
head(s)
res3 <- s %>% 
   dplyr::select(Genus, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(Genus) %>% 
   summarize(stat=mean(as.numeric(stat)))
res3

#Deframe Rank
ranks3 <- deframe(res3)
head(ranks3, 20)

#Run GSEA Analysis comparing 16S to Metagenome
fgseaRes <- fgsea(pathways=gmt.file, stats=ranks3, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #9
#Count Number of KOs in 16S Data
length(unique(s[["Genus"]])) #121
#Count Number of KOs in Metagenome Data
length(unique(a[["Genus"]])) #16

#16S (121) is Number minus the overalp and Metagenome (16) is Number minus the overlap
fit <- euler(c("16S" = 112, Metagenome = 7, "16S&Metagenome" = 9)) 

#---------------
#----Figure 2B
#---------------
pdf("16S_vs_Metagenome_TAXA_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#FF2F92","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()




#remove any NAs
s <- s[!is.na(s$Genus),]

#make sure numeric variables are numeric
s$logFC <- as.numeric(s$logFC)
s$padj <- as.numeric(s$adj.P.Val)
s <- s[!is.na(s$adj.P.Val),]

#Summarize Taxa
res4 <- s %>%
         group_by(Genus) %>% 
         summarise(KO_Number = list(unique(KO_Number)))

#Creat a List of Pathways with associated KEGG IDs
pathways.hallmark<- lapply(split(x = res4$KO_Number, f = res4$KO_Subclass_2), unlist)

#Check that the list looks good
pathways.hallmark %>% 
   head() %>% 
   lapply(head)

#Make sure it is a list
class(pathways.hallmark)

#Run GSEA Analysis Comparing 16S to Metatranscriptome
fgseaRes <- fgsea(pathways=gmt.file, stats=ranks, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #6
#Count Number of KOs in 16S Data
length(unique(s[["Genus"]])) #121
#Count Number of KOs in Metagenome Data
length(unique(x[["Genus"]])) #58

#Metatranscriptome (58) is Number minus the overalp and 16S (121) is Number minus the overlap
fit <- euler(c("16S" = 115, Metatranscriptome = 52, "16S&Metatranscriptome" = 6)) 

#---------------
#----Figure 2D
#---------------
pdf("16S_vs_Metatranscriptome_TAXA_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#FF2F92","#7A81FF"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()


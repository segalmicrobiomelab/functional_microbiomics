#Figures from the Paper
#Figure 3B

#Load Libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(eulerr)

#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Load 16S Data
s <- read.table("KEGG_PICRUST_BAL.txt", header=T, sep="\t", row.names=1)

#Load Metatranscriptome Data
x <- read.table("Metatrans_BAL_KEGG_SPT_vs_BPT 2.txt", header=T, sep="\t", row.names=1)

#Load Metagenome Data
a<-read.table("Metagenome_BAL_KEGG_SPT_vs_BPT.txt", header=T, sep="\t", row.names=1)

#//PICRUST//
#Remove NAs
data.s <- s[!is.na(s$KO_Subclass_2),]
data.s <- data.s[!is.na(data.s$logFC),]
data.s <- data.s[!is.na(data.s$adj.P.Val),]
#Calculate Median IQR and N of Padg
data.s <- setDT(data.s)[,list(Number=.N,medianfc=median(logFC), median=as.numeric(median(adj.P.Val)), iqr=as.numeric(quantile(logFC, probs=.75)),iqr1=as.numeric(quantile(logFC, probs=.25))), by=c("KO_Subclass_2")]
#Label Sequencing Type
data.s$seq <- "16S"
#convert adjusted pvalue to log10
data.s$sig <- -log10(data.s$median)
#Make sure no values become infinity
sum(is.infinite(data.s$sig))

#//Metagenome//
#Remove NAs
data.a <- a[!is.na(a$KO_Subclass_2),]
data.a <- data.a[!is.na(data.a$logFC),]
data.a <- data.a[!is.na(data.a$adj.P.Val),]
#Calculate Median IQR and N of Padg
data.a <- setDT(data.a)[,list(Number=.N,medianfc=median(logFC), median=as.numeric(median(adj.P.Val)), iqr=as.numeric(quantile(logFC, probs=.75)),iqr1=as.numeric(quantile(logFC, probs=.25))), by=c("KO_Subclass_2")]
#Label Sequencing Type
data.a$seq <- "Metagenome"
#convert adjusted pvalue to log10
data.a$sig <- -log10(data.a$median)
#Make sure no values become infinity
sum(is.infinite(data.a$sig))

#//Metatranscriptome//
#Remove NAs
data.x <- x[!is.na(x$KO_Subclass_2),]
data.x <- data.x[!is.na(data.x$log2FoldChange),]
data.x <- data.x[!is.na(data.x$padj),]
#Calculate Median IQR and N of Padg
data.x <- setDT(data.x)[,list(Number=.N,medianfc=median(log2FoldChange), median=as.numeric(median(padj)), iqr=as.numeric(quantile(log2FoldChange, probs=.75)),iqr1=as.numeric(quantile(log2FoldChange, probs=.25))), by=c("KO_Subclass_2")]
#Label Sequencing Type
data.x$seq <- "Metatranscriptome"
#convert adjusted pvalue to log10
data.x$sig <- -log10(data.x$median)
#Make sure no values become infinity
sum(is.infinite(data.x$sig))

#data.x$bb <- ifelse(data.x$sig>10, 1, 2)
#data.x$sig <- ifelse(data.x$sig>10,10,data.a$sig)

data.x <- read.table("Metatrans_BAL_KEGG_SPT_vs_BPT 3.txt", header=T, sep="\t", row.names=1)
data.x$seq <- "Metatranscriptome"
data.x$iqr1 <- NA


#keep only shared KO_Subclass
keepKOs <- data.x$KO_Subclass_2
data.s <- data.s[data.s$KO_Subclass_2 %in% keepKOs,]
keepKOs <- data.s$KO_Subclass_2
data.x <- data.x[data.x$KO_Subclass_2 %in% keepKOs,]

#merge 3 datasets
data <- rbind(data.s,data.a,data.x)

#Remove unwanted pathways
data <- data[data$KO_Subclass_2!="d pathways",]
data <- data[data$KO_Subclass_2!="Hippo signaling pathway - fly",]
data <- data[data$KO_Subclass_2!=" MAPK signaling pathway - fly",]
data <- data[data$KO_Subclass_2!="MAPK signaling pathway - yeast",]
data <- data[data$KO_Subclass_2!="Longevity regulating pathway - worm",]
data <- data[data$KO_Subclass_2!="Longevity regulating  - multiple species",]
data <- data[data$KO_Subclass_2!="Apoptosis  - multiple species",]
data <- data[data$KO_Subclass_2!=" Apoptosis - fly",]
data <- data[data$KO_Subclass_2!=" Autophagy - yeast",]
data <- data[data$KO_Subclass_2!=" MAPK   - fly",]
data <- data[data$KO_Subclass_2!=" MAPK   - plant",]


#Fix Names of Pathways
data$KO_Subclass_2 <- gsub("AGE-RAGE   in diabetic complications", "AGE-RAG (DM)",data$KO_Subclass_2)
data$KO_Subclass_2 <- gsub("Signaling  regulating pluripotency of stem cells", "Stem Cell Pluriopotency",data$KO_Subclass_2)
data$KO_Subclass_2 <- gsub("PD-L1 expression and PD-1 checkpoint  in cancer", "PD-L1",data$KO_Subclass_2)
data$KO_Subclass_2 <- gsub("-O", "O",data$KO_Subclass_2)

#change Order of dataframe to be order of 16S data by Significance
df <- data[with(data, order(seq, +sig, KO_Subclass_2)),]
#change Order of dataframe to be order of 16S data by logFC
df <- data[with(data, order(seq, +medianfc, KO_Subclass_2)),]

#Keep the Uniqe Pathways
df$KO_Subclass_2 <- factor(df$KO_Subclass_2, levels = unique(df$KO_Subclass_2))

#Set the shape of the Figures
df$shape<- as.character(as.numeric(ifelse(df$seq=="16S",21, 
	ifelse(df$seq=="Metatranscriptome",22,23))))

#Covert LOGFC to Log10
df$log10 <- df$medianfc*log10(2)
#Fix the value that is less than -4 and greater than 5
df$log10 <- ifelse(df$log10< -6, -2, df$log10)
df$log10 <- ifelse(df$log10>5, df$log10-1, df$log10)


#set variables for Colors
df$col <- ifelse(df$seq=="16S" & df$median<0.05, "A",
            ifelse(df$seq=="Metagenome" & df$median<0.05, "B",
            ifelse(df$seq=="Metatranscriptome" & df$median<0.05, "C", "D")))

#Plot Figure
pdf("All_3_Sequencing_Compared_Bubble_Chart.pdf", height = 20, width = 15)
ggplot(df, aes(x=log10, y=KO_Subclass_2, fill=col,size=Number,shape=seq)) +
	scale_shape_manual(values=c( 21,23,22), guide=FALSE) +
    geom_point(aes(fill=col),color="black",alpha=0.8)+
	#geom_point(color="black")+
    scale_size_continuous(range=c(1, 27),guide=FALSE)+
    #scale_x_reverse()+ #Flip x axis so it goes from least significant to most
    #scale_colour_gradient(low="blue", high="black")+ #Set Color for gradient
    theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(colour="black",face="bold",size=14),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
		legend.key = element_rect(colour = NA),
		legend.background = element_rect(color=NA))+
    #facet_grid(bb ~ .)+
    #scale_colour_continuous(guide = FALSE)
	scale_fill_manual(values=c( "#FF2F92","#FFD479","#7A81FF","white"), guide=FALSE) +
	#scale_color_manual(values=c( "#FF2F92","#FFD479","#7A81FF","white")) +
    #scale_x_continuous(limits=c(0,NA),expand=c(0,0),breaks=c(0,45:60), labels=c(0,45:60))+
    #annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)+
    xlab("Median Log10 Fold Change") +
    ylab("")+
    #labs( size = "Number of KOs", color="Median Log P Value" ) +
	geom_vline(xintercept=0, color="red",linetype="dashed")+
	guides(fill=FALSE)
	#theme_bw()
dev.off()




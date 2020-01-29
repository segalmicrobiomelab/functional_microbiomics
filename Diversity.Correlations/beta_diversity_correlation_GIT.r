#//////////////////////////////////////////////////////
#------------------>Histogram of Diversity
#//////////////////////////////////////////////////////
#Figures from the Paper
#Figure 5A

#Load Packages
library(grid)
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

#BetaDiversity
theme(panel.spacing = unit(2, "lines"))
beta  <- read.delim2("Beta.Diversity.txt", sep="\t")
beta$value <- as.numeric(as.character(beta$value))

#Set Order Of Figure
beta$or <-ifelse(beta$seq=="16S", 1,NA)
beta$or <-ifelse(beta$seq=="WGS",2 ,beta$or)
beta$or <-ifelse(beta$seq=="RNA",3 ,beta$or)

#Set Order of Facet
neworder <- c("ddPCR","Acetate","Propionate","Isovalerate","Butyrate")
library(plyr)  ## or dplyr (transform -> mutate)
beta <- arrange(transform(beta,
             cont=factor(cont,levels=neworder)),cont)

#Remove DDPCR
beta.nodd <- beta[beta$cont!="ddPCR",]

#--------------------
#Figure 5A
#-------------------
pdf("Beta_Diversity_Correlations_noDDPCR.pdf", height = 10, width = 25)
  ggplot(beta.nodd, aes(x=reorder(seq,+or), y=value)) + 
  geom_bar(stat="identity",color="black", width=1,
    fill=ifelse(beta.nodd$sig==1 & beta.nodd$seq=="RNA","#7A81FF",
    ifelse(beta.nodd$sig==1 & beta.nodd$seq=="16S","#FF2F92",
    ifelse(beta.nodd$sig==1 & beta.nodd$seq=="WGS","#FFD479", "white" ))))+
  facet_wrap(~cont,nrow=1)+
  xlab("")+
  ylab("PERMANOVA R2")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45))+ 
  scale_fill_manual(values=c("#FF2F92","#7A81FF","#FFD479"),guide=FALSE) +
  theme(axis.text=element_text(size=35),axis.ticks.length =unit(5,"mm"),panel.grid.major = element_blank(), panel.spacing = unit(4, "lines"), 
  panel.grid.minor = element_blank(), panel.background= element_blank(),strip.text = element_text(size=50,face="bold"), 
  axis.line = element_line(colour = "black"), strip.text.y = element_blank(), axis.text.x = element_text(),
  legend.title = element_blank(), strip.background = element_blank(), axis.title=element_text(size=30,face="bold"))
dev.off()





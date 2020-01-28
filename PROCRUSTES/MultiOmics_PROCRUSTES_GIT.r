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
library(themetagenomics)
library(grid)

##KEGG
#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),axis.title=element_text(size=20,face="bold"),
    plot.margin=unit(c(1,1,1,1),"line"),axis.text=element_text(size=15),strip.text = element_text(size=50,face="bold"),
    legend.key=element_rect(fill='white'))

#Set 16S Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/16S/")
#Load File
load(file="SCFA.DMM.RData")
#Run Vegdist for 16S Data
rld <- rld[, rld$Metatranscriptome  %in% c(1)]
rld <- rld[,colnames(rld)!="SmNV.0036.Sup.173"]
rld <- rld[,colnames(rld)!="SmNV.0038.BAL.L.171"]
rld <- rld[,colnames(rld)!="SmNV.0039.BAL.L.171"]
colnames(rld) <-gsub("171","untouched",colnames(rld))
colnames(rld) <-gsub("Sup.173","BronchSup",colnames(rld))
colnames(rld) <-gsub("Bronch.172","Bronch",colnames(rld))
colnames(rld) <-gsub("COPD0035.BronchSup","COPD.0035.BronchSup",colnames(rld))
colnames(rld) <-gsub("COPD0030.BronchSup","COPD.0030.BronchSup",colnames(rld))
colnames(rld) <-gsub("SmNV0023.BronchSup","SmNV.0023.BronchSup",colnames(rld))
colnames(rld) <-gsub("SmNV0032.BronchSup","SmNV.0032.BronchSup",colnames(rld))
colnames(rld) <-gsub("SmNV0027.BronchSup","SmNV.0027.BronchSup",colnames(rld))
colnames(rld) <-gsub("SmNV0039.BronchSup","SmNV.0039.BronchSup",colnames(rld))
rld30 <- ifelse(assay(rld)<0,0,assay(rld))

sixteens   = vegdist(t(rld30), method="bray")

CmdScale.16s <- cmdscale(sixteens, k =10)

#order by sample name
CmdScale.16s <-CmdScale.16s[order(rownames(CmdScale.16s)),]

table(rownames(CmdScale.16s)==rownames(CmdScale.trans))

#Set Metatrans Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/Metatranscriptome/")
#Load File
load(file="Metatranscriptome.RData")
#Run Vegdist for Metatrans Data
metatrans   = vegdist(t(rld30), method="bray")
CmdScale.trans <- cmdscale(metatrans, k =10)

#Set Metagenome Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/Metagenome/")
#Load File
load(file="Metagenome.RData")
#Remove data that isn't in MetaTrans
rld3 <- rld3[,rld3$SampleID!="SmNV.0039.BAL.L.untouched"]
rld3 <- rld3[,rld3$SampleID!="SmNV.0038.BAL.L.untouched"]
#Fix names so they match the MetaTrans
colnames(rld3) <- gsub("Oral", "BronchSup",colnames(rld3))
colnames(rld3) <- gsub("BALF", "BAL",colnames(rld3))
colnames(rld3) <- gsub("U", "u",colnames(rld3))
#Get data ready for Bray-Curtis Distance
rld30 <- ifelse(assay(rld3)<0,0,assay(rld3))
#Run Vegdist for Metatrans Data
metagenome   = vegdist(t(rld30), method="bray")
CmdScale.genome <- cmdscale(metagenome, k =10)

#Run the Procrustes Analysis
procrustes.results.1 <- procrustes(CmdScale.trans, CmdScale.genome)
procrustes.results.2 <- procrustes(CmdScale.16s, CmdScale.trans)
procrustes.results.3 <- procrustes(CmdScale.16s, CmdScale.genome)

#Test Significance
protest.1 <- protest(metatrans, metagenome)
protest.2 <- protest(sixteens, metatrans)
protest.3 <- protest(sixteens, metagenome)

#Metatrans versus MetaGenome
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.1$Yrot[,1],rda2=procrustes.results.1$Yrot[,2],xrda1=procrustes.results.1$X[,1],
xrda2=procrustes.results.1$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_MetaGenome_vs_MetaTrans_KEGG.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="MetaGenome"),color="#FFD479",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaTrans"),color="#7A81FF",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","MetaGenome"="#7A81FF","MetaTrans"="#FFD479","Sup"="#932CE7"),guide=guide_legend(title="Sample Type")) + 
    scale_shape_manual(values=c("MetaGenome"=16,"MetaTrans"=16), guide=guide_legend(title="Sequencing",override.aes=list(colour=c("#FFD479","#7A81FF"),size=c(3, 3))))+
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
dev.off()

#16S versus Metatrans
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.2$Yrot[,1],rda2=procrustes.results.2$Yrot[,2],xrda1=procrustes.results.2$X[,1],
xrda2=procrustes.results.2$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_16S_vs_MetaTrans_KEGG.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="16S"),color="#FF2F92",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaTrans"),color="#7A81FF",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","16S"="#7A81FF","MetaTrans"="#FF2F92","Sup"="#932CE7"),guide=guide_legend(title="Sample Type")) + 
    scale_shape_manual(values=c("16S"=16,"MetaTrans"=16), guide=guide_legend(title="Sequencing",override.aes=list(colour=c("#FF2F92","#7A81FF"),size=c(3, 3))))+
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
dev.off()

#16S versus Metagenome
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.3$Yrot[,1],rda2=procrustes.results.3$Yrot[,2],xrda1=procrustes.results.3$X[,1],
xrda2=procrustes.results.3$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_16S_vs_MetaGenome_KEGG.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="16S"),color="#FF2F92",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaGenome"),color="#FFD479",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_shape_manual(values=c("16S"=16,"MetaGenome"=16), guide=guide_legend(order=1,title="Sequencing",override.aes=list(colour=c("#FF2F92","#FFD479"),size=c(3, 3))))+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","16S"="#FFD479","MetaGenome"="#FF2F92","Sup"="#932CE7"),guide=guide_legend(order=0,title="Sample Type")) + 
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
dev.off()

#////////////////////
#------------>TAXA
#////////////////////
#Set 16S Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/16S/")
#Load File
load(file="16S.071519.RData")
#Run Vegdist for 16S Data
exp <- exp[, exp$Metatranscriptome  %in% c(1)]
rld30 <- ifelse(assay(exp)<0,0,assay(rld))
rld30 <- rld30[,colnames(rld30)!="SmNV.0036.Sup.173"]
rld30 <- rld30[,colnames(rld30)!="SmNV.0038.BAL.L.171"]
rld30 <- rld30[,colnames(rld30)!="SmNV.0039.BAL.L.171"]
colnames(rld30) <-gsub("171","untouched",colnames(rld30))
colnames(rld30) <-gsub("Sup.173","BronchSup",colnames(rld30))
colnames(rld30) <-gsub("Bronch.172","Bronch",colnames(rld30))
colnames(rld30) <-gsub("COPD0035.BronchSup","COPD.0035.BronchSup",colnames(rld30))
colnames(rld30) <-gsub("COPD0030.BronchSup","COPD.0030.BronchSup",colnames(rld30))
colnames(rld30) <-gsub("SmNV0023.BronchSup","SmNV.0023.BronchSup",colnames(rld30))
colnames(rld30) <-gsub("SmNV0032.BronchSup","SmNV.0032.BronchSup",colnames(rld30))
colnames(rld30) <-gsub("SmNV0027.BronchSup","SmNV.0027.BronchSup",colnames(rld30))
colnames(rld30) <-gsub("SmNV0039.BronchSup","SmNV.0039.BronchSup",colnames(rld30))
sixteens   = vegdist(t(rld30), method="bray")

CmdScale.16s <- cmdscale(sixteens, k =10)

#order by sample name
CmdScale.16s <-CmdScale.16s[order(rownames(CmdScale.16s)),]

#table(rownames(CmdScale.16s)==rownames(CmdScale.trans))

#Set Metatrans Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/Metatranscriptome/")
#Load File
load(file="Metatranscriptome.RData")
#Run Vegdist for Metatrans Data
metatrans   = vegdist(t(rld20), method="bray")
CmdScale.trans <- cmdscale(metatrans, k =10)

#Set Metagenome Working Directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/Metagenome/")
#Load File
load(file="Metagenome.RData")
#Remove data that isn't in MetaTrans
rld2 <- rld2[,rld2$SampleID!="SmNV.0039.BAL.L.untouched"]
rld2 <- rld2[,rld2$SampleID!="SmNV.0038.BAL.L.untouched"]
#Fix names so they match the MetaTrans
colnames(rld2) <- gsub("Oral", "BronchSup",colnames(rld2))
colnames(rld2) <- gsub("BALF", "BAL",colnames(rld2))
colnames(rld2) <- gsub("U", "u",colnames(rld2))
#Get data ready for Bray-Curtis Distance
rld20 <- ifelse(assay(rld2)<0,0,assay(rld2))
#Run Vegdist for Metatrans Data
metagenome   = vegdist(t(rld20), method="bray")
CmdScale.genome <- cmdscale(metagenome, k =10)

#Run the Procrustes Analysis
procrustes.results.1 <- procrustes(CmdScale.trans, CmdScale.genome)
procrustes.results.2 <- procrustes(CmdScale.16s, CmdScale.trans)
procrustes.results.3 <- procrustes(CmdScale.16s, CmdScale.genome)

#Test Significance
protest.1 <- protest(metatrans, metagenome)
protest.2 <- protest(sixteens, metatrans)
protest.3 <- protest(sixteens, metagenome)

#Metatrans versus MetaGenome
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.1$Yrot[,1],rda2=procrustes.results.1$Yrot[,2],xrda1=procrustes.results.1$X[,1],
xrda2=procrustes.results.1$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_MetaGenome_vs_MetaTrans_TAXA.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="MetaGenome"),color="#FFD479",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaTrans"),color="#7A81FF",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","MetaGenome"="#7A81FF","MetaTrans"="#FFD479","Sup"="#932CE7"),guide=guide_legend(title="Sample Type")) + 
    scale_shape_manual(values=c("MetaGenome"=16,"MetaTrans"=16), guide=guide_legend(title="Sequencing",override.aes=list(colour=c("#FFD479","#7A81FF"),size=c(3, 3))))+
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
    dev.off()




#16S versus Metatrans
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.2$Yrot[,1],rda2=procrustes.results.2$Yrot[,2],xrda1=procrustes.results.2$X[,1],
xrda2=procrustes.results.2$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_16S_vs_MetaTrans_TAXA.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="16S"),color="#FF2F92",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaTrans"),color="#7A81FF",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","16S"="#7A81FF","MetaTrans"="#FF2F92","Sup"="#932CE7"),guide=guide_legend(title="Sample Type")) + 
    scale_shape_manual(values=c("16S"=16,"MetaTrans"=16), guide=guide_legend(title="Sequencing",override.aes=list(colour=c("#FF2F92","#7A81FF"),size=c(3, 3))))+
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
dev.off()

#16S versus Metagenome
#Get Data Ready for Plotting
ctest <- data.frame(rda1=procrustes.results.3$Yrot[,1],rda2=procrustes.results.3$Yrot[,2],xrda1=procrustes.results.3$X[,1],
xrda2=procrustes.results.3$X[,2])

ctest <- merge(ctest,colData(rld2),by=0)
ctest <- as.data.frame(ctest)

pdf("Procrustes_16S_vs_MetaGenome_TAXA.pdf", height = 10, width = 10)
    ggplot(ctest) +
    #geom_point(aes(x=xrda1, y=xrda2,fill="MetaGenome",color="MetaGenome"),size=3) +
    #geom_point(aes(x=rda1, y=rda2,fill="MetaTrans",color="MetaTrans"),size=3) +
    geom_point(aes(x=xrda1, y=xrda2,shape="16S"),color="#FF2F92",size=5,alpha=0.8) +
    geom_point(aes(x=rda1, y=rda2,shape="MetaGenome"),color="#FFD479",size=5,alpha=0.8) +
    #geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),arrow=arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(xend=rda1,yend=rda2,x=xrda1,y=xrda2,color=Sample_Type_DMM_Class_ReSeq),alpha=0.6)+
    xlab("Axis 1")+
    ylab("Axis 2")+
    scale_shape_manual(values=c("16S"=16,"MetaGenome"=16), guide=guide_legend(order=1,title="Sequencing",override.aes=list(colour=c("#FF2F92","#FFD479"),size=c(3, 3))))+
    scale_color_manual(values=c("BAL.BPT"="#296218","BAL.SPT"="#EA3323","BKG"="#000000","16S"="#FFD479","MetaGenome"="#FF2F92","Sup"="#932CE7"),guide=guide_legend(order=0,title="Sample Type")) + 
    #scale_shape_manual(labels=c("MetaGenome","MetaTrans"),values=c(19,17))+
    theme
dev.off()


save.image(file="Proscrustes.RData")

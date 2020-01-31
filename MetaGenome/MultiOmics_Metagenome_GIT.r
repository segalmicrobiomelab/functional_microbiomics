#-----------------------------------------
#---------------------->Metagenome Analysis
#-----------------------------------------

#Figures from the Paper
#Figure 1G
#SFigure 1A
#Figure 2C
#Figure 5C
#SFigure 5A
#SFigure 8B
#SFigure 9B
#SFigure 8A


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
library(vegan)
library(ggpmisc)
library(dplyr)



#Set Theme for Figures
theme<-	theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.05

#Load Meta Data
coldata <- read.delim2("Map.Metatrancriptome.COPD.SmNV.b2.txt", sep="\t")

#Order Meta Data by SampleId
coldata <- coldata[order(coldata$SampleID),]

#load Count Data
mycounts <-read.delim2("Metagenome_KEGG_ANNOTATED.txt", sep="\t", row.names=1)

#Order Count Data by SampleID
mycounts <-mycounts[, order(colnames(mycounts))]

#Confirm Sample IDs match for Count and Meta Data
table(colnames(mycounts)==as.character(coldata$SampleID))

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0

#Copy of Count Table
mycounts3 <- mycounts

#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts3, function(x) as.numeric(as.character(x))),
                   check.names=F, row.names = rownames(mycounts3))


#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

#Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata,
                              design= ~ Sample_Type_DMM_Class_ReSeq)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds), 1, gm_mean)

#Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- estimateSizeFactors(dds)

#Transforming data - Option #1 is regularized-logarithm transformation, or rlog for short. 
#For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation. 
#For genes with lower counts, however, the values are shrunken towards the genes' averages across all sample
rld <- rlog(dds, fitType="local",blind=FALSE)

#Transforming data - Option #2 is variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds)

#keep only data that is Annotated to a KEGG
dds1 <- dds[!grepl("UNGROUPED|UNMAPPED",rownames(assay(dds))),]
rld1 <- rld[!grepl("UNGROUPED|UNMAPPED",rownames(assay(rld))),]
vsd1 <- vsd[!grepl("UNGROUPED|UNMAPPED",rownames(assay(vsd))),]

#keep only data with Taxa Annotation
dds2 <- dds1[grep("g_|unclassified",rownames(assay(dds1))),]
dds2 <- dds2[!grepl("unclassified",rownames(assay(dds2))),]
rld2 <- rld1[grep("g_|unclassified",rownames(assay(rld1))),]
rld2 <- rld2[!grepl("unclassified",rownames(assay(rld2))),]
vsd2 <- vsd1[grep("g_|unclassified",rownames(assay(vsd1))),]
vsd2 <- vsd2[!grepl("unclassified",rownames(assay(vsd2))),]

#keep only data for KEGG
dds3 <- dds1[!grepl("g_|unclassified",rownames(assay(dds1))),]
rld3 <- rld1[!grepl("g_|unclassified",rownames(assay(rld1))),]
vsd3 <- vsd1[!grepl("g_|unclassified",rownames(assay(vsd1))),]

#===========================================================================================================================================================================
#===========================================================================================================================================================================
#////////////////// KEGG /////////////////////
#===========================================================================================================================================================================
#===========================================================================================================================================================================

#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================
#Create Distance Matrix
#vegdist   = vegdist(t(assay(rld3)), method="bray")
rld30 <- ifelse(assay(rld3)<0,0,assay(rld3))
vegdist   = vegdist(t(rld30), method="bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(rld3), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Sample_Type_DMM_Class_ReSeq,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_DMM_Class_ReSeq",suffixes=c("",".centroid"))

#-----------
#Figure 1G
#-----------
pdf("Metagenome_KEGG_Sample_Type_DMM_Class_BRAY.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Sample_Type_DMM_Class_ReSeq)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Sample_Type_DMM_Class_ReSeq), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Sample_Type_DMM_Class_ReSeq))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL.BPT", "BAL.SPT", "BKG", "Sup")), size=8) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Check Statistics
adonis(vegdist ~ vsd3$Sample_Type_DMM_Class_ReSeq)
                                 Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
vsd3$Sample_Type_DMM_Class_ReSeq  3   0.52151 0.173837  5.0772 0.30939  0.001
Residuals                        34   1.16412 0.034239         0.69061
Total                            37   1.68563                  1.00000



#Re Do Figure with Clustering based on MetaGenome
#Merge PC Data with MetaData
require(data.table)
rownames(coldata) <- coldata$SampleID
newResults <- merge(x = CmdScale, y = coldata, by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Sample_Type_DMM_Class_ReSeq_Meagenome,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_DMM_Class_ReSeq_Meagenome",suffixes=c("",".centroid"))

#-----------
#SFigure 9B
#-----------
pdf("Metagenome_KEGG_Sample_Type_Metagenome_DMM_Class_BRAY.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Sample_Type_DMM_Class_ReSeq_Meagenome)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Sample_Type_DMM_Class_ReSeq_Meagenome), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Sample_Type_DMM_Class_ReSeq_Meagenome))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL.BPT", "BAL.SPT", "BKG", "UA")), size=8) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#=========================================================
/////////////////////Alpha Diversity///////////////////////
#=========================================================

#Calcultes Shannon Diversity
dds3$Shannon = diversity(rld30, index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
shannon = as.data.frame(colData(dds3))
#Remove any zero values
shannon[shannon==0] <- NA

#Set Order Of Figure
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BKG", 1,NA)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BAL.BPT",2 ,shannon$or)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BAL.SPT",3 ,shannon$or)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="Sup",4 ,shannon$or)

#-----------
#SFigure 1A
#-----------
pdf("Metagenome_KEGG_Sample_Type_DMM_Class_SHANNON.pdf", height = 7, width = 5)
    ggplot(shannon, aes(x= reorder(Sample_Type_DMM_Class_ReSeq, +or), y=Shannon, color=Sample_Type_DMM_Class_ReSeq)) + 
    stat_boxplot(geom ='errorbar', width=0.1)+
    geom_boxplot(aes(color=Sample_Type_DMM_Class_ReSeq),outlier.shape = NA, width=0.5)+
    #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
    geom_jitter(aes(fill=Sample_Type_DMM_Class_ReSeq), position=position_jitter(0.2))+
    scale_fill_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) +
    scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) +
    scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
    #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("Shannon Diversity") + 
    xlab("Sample Type")+
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme
dev.off()

#Check Statistics
kruskal.test(Shannon ~ Sample_Type_DMM_Class_ReSeq, data = shannon)
Kruskal-Wallis chi-squared = 1.6, df = 3, p-value = 0.6594


#=========================================================
//////////InterGroup Beta Diversity///////////////////////
#=========================================================
#######################   DON'T CHANGE THIS FUNCTION: STARTS HERE     ##################################################
intergroup.distances <- function(sample.data, Variable.Intergroup, distance.matrix, extraVar.addtotable, filename){
#melt the distance matrix into columns with each sample site and distance
b <- melt(as.matrix(distance.matrix))

#Then need to remove self distances and duplicated distances
p    <- t(apply(b[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])

p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))

#establish new data frame that removes those values 
b.df   <- b[-c(rmv1,rmv2),] 

##Now we need to ADD variable rows with location 
#set up new data frame
new.df <- b.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new.df[] <- lapply(b.df, function(x) Variable.Intergroup[match(x, rownames(sample.data))])

#create two lists of the group variable 
topo.var1 <- new.df[,1]
topo.var2 <-new.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b.var.df <- cbind(b.df, topo.var1, topo.var2)

##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b.var.df <- b.var.df
#set row names to re-zero
rownames(btw.b.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b.var.df)) {
if (btw.b.var.df$topo.var1[i] == btw.b.var.df$topo.var2[i]) {
toremove <- append(toremove, i)

} 
}

#remove indexes we selected
btw.b.var.df <- btw.b.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw = paste(btw.b.var.df$topo.var1, "to", btw.b.var.df$topo.var2)
new.cat.btw.df <- data.frame(btw.b.var.df$topo.var1, btw.b.var.df$topo.var2)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw.df, 1, sort))
unique.new.cat.btw <- unique(new.cat.btw.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw) <- NULL
rownames(unique.new.cat.btw) <- NULL
unique.new.cat.btw <- paste(unique.new.cat.btw[,1], "to", unique.new.cat.btw[,2])

#create new data frame 
clean.btw.b.var.df <- btw.b.var.df

#reset row names
rownames(clean.btw.b.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b.var.df)){
if (paste(clean.btw.b.var.df$topo.var2[i], "to", clean.btw.b.var.df$topo.var1[i]) %in% unique.new.cat.btw) {
clean.btw.b.var.df$topo.var1[i] <- btw.b.var.df$topo.var2[i]
clean.btw.b.var.df$topo.var2[i] <- btw.b.var.df$topo.var1[i]
}
}

#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw.clean = paste(clean.btw.b.var.df$topo.var1, "to", clean.btw.b.var.df$topo.var2)

#confirm permutations 
unique(new.cat.btw.clean)

# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
comb.inter.data.braypart <- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value))

##Now we need to ADD variable rows with location 
var1 <- comb.inter.data.braypart$Sample1
var2 <- comb.inter.data.braypart$Sample2

x <- match(var1, rownames(sample.data))
y <- extraVar.addtotable[x]
#y now has the extra var for var 1 


x2 <- match(var2, rownames(sample.data))
y2 <- extraVar.addtotable[x2]
#y2 has the extra var for var 2


##Add the category 
b <- Variable.Intergroup[x]
#b now has the grouping var for var 1 

b2 <- Variable.Intergroup[x2]
#b2 has the grouping var var 2


# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs AND the subject IDs for comparison later 
allEncompassing.df <<- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value), Subject1 = y, SampleType1 = b, Subject2 = y2, SampleType2 = b2)
}
#######################   DON'T CHANGE THIS FUNCTION: ENDS HERE     ##################################################


#-------->//Sup to BAL// vs SCFA
#Create the tables for the function
Sup.BAL <- rld3[, rld3$Sample_Type_Simple!= "BKG"]
Sup.BAL.mat <- ifelse(assay(Sup.BAL)<0,0,assay(Sup.BAL))
Bray.dist   = vegdist(t(Sup.BAL.mat), method="bray")


######## Here, you set your parameters. In this case we use:   Bray  ###################
{
#set the sample data 
sample.data = colData(Sup.BAL)
#set the variable you want to do groupings by
Variable.Intergroup = colData(Sup.BAL)$Sample_Type_DMM_Class_ReSeq
#set the distance.matrix
distance.matrix = Bray.dist
##Select if you would like to add any other variables to the final table for clarifiation and processing. If none, then just repeat the variable.intergroup
#in this care s we need to add the Subject_ID, in this case: ID_Sample_Type_Subject_Type_Simple 
extraVar.addtotable = colData(Sup.BAL)$Sample_Type_DMM_Class_ReSeq
#set the file name you'd like
filename = "Unique.ID.Bronch.Cohort.Bray.txt"
allEncompassing.df = NULL


#run it for all distances between samples
intergroup.distances(sample.data, Variable.Intergroup, distance.matrix, extraVar.addtotable, filename)


#Remove any BAL to BAL comparissons
allEncompassing.df <- allEncompassing.df[allEncompassing.df$Category!="BAL.SPT to BAL.BPT",]
allEncompassing.df <- allEncompassing.df[allEncompassing.df$Category!="BAL.BPT to BAL.SPT",]

#Create a Variable that is BAL to BKG Comparisson
allEncompassing.df$comparison <- ifelse(allEncompassing.df$Subject1=="Sup",
                                paste(allEncompassing.df$Sample2,"to",allEncompassing.df$Sample1), paste(allEncompassing.df$Sample1,"to",allEncompassing.df$Sample2))
#Remove Duplicates
allEncompassing.df = allEncompassing.df[!duplicated(allEncompassing.df$comparison),]

#Create a variable for SampleID
allEncompassing.df$SampleID <- ifelse(allEncompassing.df$Subject1=="Sup",
                                as.character(allEncompassing.df$Sample2), as.character(allEncompassing.df$Sample1))

#Keep only the Sample ID and the distance
allEncompassing.df <- allEncompassing.df[,names(allEncompassing.df) %in% c("SampleID","Distance")]

#Summarize the Mean of the data per sampleID
data <- allEncompassing.df %>%
  group_by(SampleID) %>%
  summarize(mean_size = mean(Distance, na.rm = TRUE))


BAL.Sample <- Sup.BAL[, Sup.BAL$Sample_Type_Simple %in% "BAL"]
data.sup <- data.frame(data,Sample_Type=colData(BAL.Sample)$Sample_Type_DMM_Class_ReSeq,
                    Acetate    =as.numeric(as.character(colData(BAL.Sample)$Acetate)),
                    Propionate =as.numeric(as.character(colData(BAL.Sample)$Propionate)),
                    Butyrate   =as.numeric(as.character(colData(BAL.Sample)$Butyrate)),
                    Isovalerate=as.numeric(as.character(colData(BAL.Sample)$Isovalerate)),
                    Valerate   =as.numeric(as.character(colData(BAL.Sample)$Valerate)),
                    Hexanoate  =as.numeric(as.character(colData(BAL.Sample)$Hexanoate)),
                    Octanoate  =as.numeric(as.character(colData(BAL.Sample)$Octanoate))
                    )


#Write Table For Figure 1I
write.table(data.sup, file = "Metagenome.Sup.BAL.Bray.txt", sep = "\t", row.names = FALSE)


#===========================================================================================================================================================================
#===========================================================================================================================================================================
#//////////////////BAL TAXA/////////////////////
#===========================================================================================================================================================================
#===========================================================================================================================================================================
#Subset BAL of TAXA Data and drop unwanted levels
dds4 <- dds2[, dds2$Sample_Type_Simple %in% "BAL"]
vsd4 <- vsd2[, vsd2$Sample_Type_Simple %in% "BAL"]
rld4 <- rld2[, rld2$Sample_Type_Simple %in% "BAL"]

#Drop Unencessary Levels to Run Differential Expression
dds4$Sample_Type_DMM_Class_ReSeq <- droplevels(dds4$Sample_Type_DMM_Class_ReSeq)

#Differential Analysis of BAL TAXA
#Set the baseline level --> positive is upregulated in BAL.SPT ; Negative is down regulated
dds4$Sample_Type_DMM_Class_ReSeq <- relevel(dds4$Sample_Type_DMM_Class_ReSeq, ref ="BAL.BPT")

#Run the differential Analysis
dds4<- DESeq(dds4, test="Wald", fitType="local")

#Create A table of differential expression analysis
res4 <- results(dds4, cooksCutoff=FALSE)


# Reorder Results based on FDR
res4 = res4[order(res4$padj, na.last = NA), ]

#=========================================================
//////////////////////TABLES/////////////////////////////
#=========================================================
#Convert Resuts table into a data.frame
res4 <- as.data.frame(res4)
#Get Genus and Species from name
res4$gs<- substring(rownames(res4), 8)
#Get KO from name
res4$ko<- substr(rownames(res4), 1,6)
#Get Species from gs
res4$species <- lapply(strsplit(as.character(res4$gs), "\\."), "[", 2)
#Get Genus from gs
res4$genus <- lapply(strsplit(as.character(res4$gs), "\\."), "[", 1)

######get abundance data 
#decide what IDs to save 
ko.to.save <-as.character(rownames(res4))

#from main table we should get the mean across the row of the table
BAL.ko.table.df <- data.frame(assay(dds4))
BAL.ko.table.df.meanRA <- rowMeans(BAL.ko.table.df)

#need to subset AND reorder just the IDs that we have 
BAL.ko.table.df.meanRA.save <- BAL.ko.table.df.meanRA[ko.to.save]

#add the abundnace data for the res dataframe
res4$abundance <- BAL.ko.table.df.meanRA.save

#Set Names of Results Table
res4 <- setNames(cbind(rownames(res4), res4, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "gs", "ko","species","genus","abundance")) 

#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================
# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res4$sig <- -log10(res4$adj.P.Val)
#Identify any infinite values created by log10
sum(is.infinite(res4$sig))
#Set a celling for any infinite values
res4[is.infinite(res4$sig),"sig"] <- 350

#Set the color density of the dots
cols <- densCols(res4$logFC, res4$sig)
cols[res4$pvalue ==0] <- "purple"
cols[res4$logFC > 2 & res4$adj.P.Val < alpha ] <- "red"

#-----------
#Figure 2C
#-----------
pdf(file="TAXA_Metagenome_BAL_SPT_vs_BPT_Volcano_Plot_Small_FDR_0.05.pdf", width=5, height=5)
    ggplot(res4, aes(x = logFC, y = sig,label=Gene.symbol)) +
    geom_point(color=cols, size = ifelse(res4$adj.P.Val < alpha, res4$abundance, 2), alpha=0.5) + #Chose Colors and size for dots
    geom_text_repel(aes(label=ifelse(res4$logFC>0 & res4$adj.P.Val < 0.05 , as.character(res4$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
    theme(legend.position = "none") +
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
    xlab("Effect size: log2(fold-change)") +
    ylab("-log10(adjusted p-value)") + 
    ylim(0,20)+
    theme
dev.off() 


#===========================================================================================================================================================================
#===========================================================================================================================================================================
#//////////////////BAL KEGG/////////////////////
#===========================================================================================================================================================================
#===========================================================================================================================================================================
#Subset BAL 
dds5 <- dds3[, dds3$Sample_Type_Simple %in% "BAL"]

#Drop Unencessary Levels to Run Differential Expression
dds5$Sample_Type_DMM_Class_ReSeq <- droplevels(dds5$Sample_Type_DMM_Class_ReSeq)

#Differential Analysis of BAL TAXA
#Set the baseline level --> positive is upregulated in BAL.SPT ; Negative is down regulated
dds5$Sample_Type_DMM_Class_ReSeq <- relevel(dds5$Sample_Type_DMM_Class_ReSeq, ref ="BAL.BPT")

#Run the differential Analysis
dds5<- DESeq(dds5, test="Wald", fitType="local")
#dds5<- DESeq(dds5,fitType="mean")

#Create A table of differential expression analysis
res4 <- results(dds5, cooksCutoff=FALSE)
res4 <- results(dds5)

# Reorder Results based on FDR
res4 = res4[order(res4$padj, na.last = NA), ]


#=========================================================
/////////////////////////Beta Diversity ///////////////////////
#=========================================================
#Create Distance Matrix
rld50 <- ifelse(assay(rld5)<0,0,assay(rld5))
vegdist   = vegdist(t(rld50), method="bray")

#-----------------------------------------------------
##Correlation Analysis with SCFA and Diversity Metrics
#-----------------------------------------------------
#Statistics with Acetate
data.adonis <- data.frame(colData(rld5))
adonis(vegdist ~ as.numeric(as.character(Acetate)), data.adonis)
                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Acetate))  1   0.62421 0.62421  7.9937 0.31983  0.001
Residuals                         17   1.32749 0.07809         0.68017
Total                             18   1.95170                 1.00000

as.numeric(as.character(Acetate)) ***
Residuals
Total

#Alpha
Shannon = diversity(rld50, index = "shannon", MARGIN = 2, base = exp(1))
#Convert to data frame for ggplot
shannon = data.frame(colData(rld5))
#Adonis
adonis(Shannon ~ as.numeric(as.character(Acetate)), shannon, na.rm=TRUE)
                                  Df SumsOfSqs    MeanSqs  F.Model      R2
as.numeric(as.character(Acetate))  1 0.0000150 0.00001496 0.039025 0.00229
Residuals                         17 0.0065168 0.00038334          0.99771
Total                             18 0.0065318                     1.00000
                                  Pr(>F)
as.numeric(as.character(Acetate))  0.861

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Acetate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

S = 1466, p-value = 0.2345
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.2859649

#Statistics with Propionate
adonis(vegdist ~ as.numeric(as.character(Propionate)), data.adonis)
                                     Df SumsOfSqs MeanSqs F.Model      R2
as.numeric(as.character(Propionate))  1   0.51263 0.51263  6.0558 0.26266
Residuals                            17   1.43907 0.08465         0.73734
Total                                18   1.95170                 1.00000
                                     Pr(>F)
as.numeric(as.character(Propionate))  0.001 ***
Residuals
Total


#Alpha Adonis
adonis(Shannon ~ as.numeric(as.character(Propionate)), shannon, na.rm=TRUE)
                                     Df SumsOfSqs    MeanSqs  F.Model      R2
as.numeric(as.character(Propionate))  1 0.0000165 0.00001652 0.043095 0.00253
Residuals                            17 0.0065153 0.00038325          0.99747
Total                                18 0.0065318                     1.00000
                                     Pr(>F)
as.numeric(as.character(Propionate))  0.816
Residuals
Total

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Propionate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
         
S = 1396, p-value = 0.3538
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.2245614

#Statistics with Isovalerate
adonis(vegdist ~ as.numeric(as.character(Isovalerate)), data.adonis)
                                      Df SumsOfSqs MeanSqs F.Model      R2
as.numeric(as.character(Isovalerate))  1   0.61372 0.61372  7.7977 0.31445
Residuals                             17   1.33798 0.07870         0.68555
Total                                 18   1.95170                 1.00000
                                      Pr(>F)
as.numeric(as.character(Isovalerate))  0.001 ***
Residuals
Total

#Alpha Adonis
adonis(Shannon ~ as.numeric(as.character(Isovalerate)), shannon, na.rm=TRUE)
                                      Df SumsOfSqs    MeanSqs  F.Model      R2
as.numeric(as.character(Isovalerate))  1 0.0000238 0.00002381 0.062187 0.00364
Residuals                             17 0.0065080 0.00038282          0.99636
Total                                 18 0.0065318                     1.00000
                                      Pr(>F)
as.numeric(as.character(Isovalerate))  0.795
Residuals
Total

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Isovalerate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
         
S = 1240, p-value = 0.7209
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.0877193



#Statistics with Butyrate
adonis(vegdist ~ as.numeric(as.character(Butyrate)), data.adonis)
                                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Butyrate))  1   0.27629 0.276288  2.8034 0.14156  0.056
Residuals                          17   1.67541 0.098554         0.85844
Total                              18   1.95170                  1.00000

#Alpha Adonis
adonis(Shannon ~ as.numeric(as.character(Butyrate)), shannon, na.rm=TRUE)
                                   Df SumsOfSqs    MeanSqs F.Model      R2
as.numeric(as.character(Butyrate))  1 0.0006709 0.00067094  1.9461 0.10272
Residuals                          17 0.0058609 0.00034476         0.89728
Total                              18 0.0065318                    1.00000
                                   Pr(>F)
as.numeric(as.character(Butyrate))  0.148
Residuals
Total

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Butyrate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
         
S = 1460, p-value = 0.2435
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.2807018


#=========================================================
//////////////////////TABLES/////////////////////////////
#=========================================================
#Convert Resuts table with pathway into a data.frame
res4 <- as.data.frame(res4)
#Create Column of Rownames
res4$KO_Number <- rownames(res4)
#Keep only the KO Number
res4$KO_Number <- substr(res4$KO_Number, 0, 6)
#Asign Hierarchy to KO Number
res4 <- assign_hierarchy(count_data=res4, keep_unknowns=TRUE, identifier="KO_Number")
#Remove Number from Subclass
res4$KO_Subclass_2 <- gsub('[[:digit:]]+', '', res4$KO_Subclass_2 )
#Remove Path Number from Subclass
res4$KO_Subclass_2 <- gsub("\\s*\\[[^\\)]+\\]","",res4$KO_Subclass_2)


#decide what otu to save 
ko.to.save <-as.character(rownames(res4))

#from relative table we should get the mean across the row of the otu table
BAL.ko.table.df <- data.frame(assay(dds5))
BAL.ko.table.df.meanRA <- rowMeans(BAL.ko.table.df)

#need to subset AND reorder just the otus that we have 
BAL.ko.table.df.meanRA.save <- BAL.ko.table.df.meanRA[ko.to.save]

#add the abundnace data for the res dataframe
res4$abundance <- BAL.ko.table.df.meanRA.save

#Set Names of Results Table
res4 <- setNames(cbind(rownames(res4), res4, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "KO_Number","KO_Class", "KO_Subclass_1","KO_Subclass_2","abundance")) 
res4$path <- paste(res4$KO_Number,res4$KO_Class,res4$KO_Subclass_1,res4$KO_Subclass_2, sep="_")


#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================
#Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res4$sig <- -log10(res4$adj.P.Val)
sum(is.infinite(res4$sig))
res4[is.infinite(res4$sig),"sig"] <- 350

#Remove NAs
res4 <- res4[!is.na(res4$KO_Subclass_2),]


#Set Density Colors
cols <- densCols(res4$logFC, res4$sig)
cols[res4$pvalue ==0] <- "purple"
cols[res4$logFC > 2 & res4$adj.P.Val < alpha ] <- "red"
res4$pch <- 19
res4$pch[res4$pvalue ==0] <- 6

#---------------------
#SFigure 5A
#---------------------
pdf(file="KEGG_Metagenome_BAL_SPT_vs_BPT_Volcano_Plot_Small_FDR_0.05.pdf", width=5, height=5)
    ggplot(res4, aes(x = logFC, y = sig,label=KO_Subclass_2)) +
    geom_point(color=cols, size = ifelse(res4$adj.P.Val < alpha, 10^-1*res4$abundance, 2), alpha=0.5) + #Chose Colors and size for dots
    geom_text_repel(aes(label=ifelse(res4$logFC>7.3 & res4$adj.P.Val < alpha , as.character(res4$path),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.5) +
    theme(legend.position = "none") +
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
    xlab("Effect size: log2(fold-change)") +
    ylab("-log10(adjusted p-value)") + 
    ylim(0,20)+
    theme
dev.off() 


#------------------------------------------------
#---->Using SCFA analysis with MetaGenome data
#------------------------------------------------
#Remove any data that is not annotated
d2 <- d1[!grepl("UNGROUPED|UNMAPPED",rownames(d1)),]

#keep only data for KEGG
kegg <- d2[!grepl("g_|unclassified",rownames(d2)),]
#Create Column for KO
kegg$ko <- substr(rownames(kegg), 1,6)

#Save KEGG Table as another object
kegg.df <- kegg

#Drop Reddundant Columns
drops <- c("ko")
kegg2 <- kegg.df[ , !(names(kegg.df) %in% drops)]

#Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds.keggscfa <- DESeqDataSetFromMatrix(countData = kegg2,
                              colData = coldata,
                              design= ~ Sample_Type_DMM_Class_ReSeq_Meagenome)
                              
out <-colSums(counts(dds.keggscfa))
out #count how many samples have 0 counts = 6; run the following code 6 times

#Remove all the 0s
idx <- which.max(colSums(counts(dds.keggscfa) == 0))
dds2 <- dds.keggscfa[ , -idx]

idx <- which.max(colSums(counts(dds2) == 0))
dds2 <- dds2[ , -idx]

idx <- which.max(colSums(counts(dds2) == 0))
dds2 <- dds2[ , -idx]

idx <- which.max(colSums(counts(dds2) == 0))
dds2 <- dds2[ , -idx]

idx <- which.max(colSums(counts(dds2) == 0))
dds2 <- dds2[ , -idx]

idx <- which.max(colSums(counts(dds2) == 0))
dds2 <- dds2[ , -idx]
                              
#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)

#Estimate Factors of DESeq Object
dds.keggscfa <- estimateSizeFactors(dds2, geoMeans = geoMeans)

#=========================================================
//////////////DOT PLOTS OF KOs///////////////////////
#=========================================================
#Convert Normalized Assay Data into a Data Table
as <- data.table(melt(assay(dds.keggscfa)))

#Create a Variable name to match with MetaData
as$SampleID <- as$Var2

#Create a Data Frame of MetaData
assay <- as.data.frame(colData(dds2))
assay$SampleID <- rownames(assay)

#Merge Assay Data and Meta Data
as2 <- merge(as,assay, all=TRUE)

# convert Genus to a character vector from a factor
as2$KEGG <- as.character(as2$Var1)

#Create Data Frame copy
as3 <- as2

#Extract KEGGs for Acetate, Propionate and Butyrate
keepKOSC <- c("K01738","K00925","K01034")
as2 <- as3[as3$KEGG %in% keepKOSC,]

#Test Overall Significance
kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2) 
Kruskal-Wallis chi-squared = 0.74144, df = 3, p-value = 0.8634

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq_Meagenome %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2.BPT.SPT)
W = 127, p-value = 0.5135

#Replace 0 with 0.1
as2$value <- ifelse(as2$value==0,0.1,as2$value)
#Change BKG label so that it comes first in order
as2$Sample_Type_DMM_Class_ReSeq_Meagenome <- as.character(as2$Sample_Type_DMM_Class_ReSeq_Meagenome)
as2$Sample_Type_DMM_Class_ReSeq_Meagenome <- ifelse(as2$Sample_Type_DMM_Class_ReSeq_Meagenome=="BKG","aBKG",as2$Sample_Type_DMM_Class_ReSeq_Meagenome)
#Change Order Of Graph
as2$order <-ifelse(as2$KEGG=="K01738",1,
            ifelse(as2$KEGG=="K00925",2,3))

#---------------------
#Figure 5C
#---------------------
pdf("Metagenome_KEGG_available_BAL_SPT_vs_BPT_BOX_Plot.pdf.pdf", height = 7, width = 10)
    ggplot(as2, aes(x=reorder(KEGG,+order), y=value, color=Sample_Type_DMM_Class_ReSeq_Meagenome)) +
    geom_boxplot(aes(color=Sample_Type_DMM_Class_ReSeq_Meagenome),outlier.shape = NA, width=0.7)+
    geom_point(aes(x=reorder(KEGG,+order), y=value),pch=21, fill="black",
             position=position_jitterdodge(jitter.width=0.1, dodge.width=0.7), size=3, alpha=1)+  
    scale_fill_manual(values=c("black"="black")) +
    scale_y_continuous(name="KEGG Expression",trans="log10", breaks = trans_breaks('log10', function(x) 10^x),limits=c(10^-1,10^5), labels = trans_format('log10', math_format(10^.x))) +
    xlab("") +
    scale_color_manual(label=c("BKG","BAL.BPT","BAL.SPT","Sup"),
        values=c("#000001","#296218","#EA3323","#932CE7"),
        guide=guide_legend(order=0,title="Sample Type")) + 
        guides(fill=FALSE)+
    scale_x_discrete(labels = c('K01738\n (Acetate)','K00925\n (Propionate)','K01034\n (Butyrate)'))+ 
    theme(axis.text=element_text(size=25),panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background= element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(), 
    strip.background = element_blank(), axis.title=element_text(size=30,face="bold"))
dev.off()

#Check Significance for Each SCFA individually
#Extract Acetate Data
keepKOSC <- c("K01738")
as2 <- as3[as3$KEGG %in% keepKOSC,]
            
kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2) 
Kruskal-Wallis chi-squared = 11.614, df = 3, p-value = 0.008831

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq_Meagenome %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2.BPT.SPT)
W = 0, p-value = 0.0636

#Extract Propionate Data
keepKOSC <- c("K00925")
as2 <- as3[as3$KEGG %in% keepKOSC,]
            
kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2) 
Kruskal-Wallis chi-squared = 13.785, df = 3, p-value = 0.003213

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq_Meagenome %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2.BPT.SPT)
W = 0, p-value = 0.0636

#Extract Butyrate Data
keepKOSC <- c("K01034")
as2 <- as3[as3$KEGG %in% keepKOSC,]
            
kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2) 
Kruskal-Wallis chi-squared = 5.3359, df = 3, p-value = 0.1488

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq_Meagenome %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq_Meagenome, data = as2.BPT.SPT)
W = 4.5, p-value = NA
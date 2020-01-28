#-----------------------------------------
#----------------------------->16S Analysis
#-----------------------------------------

#Figures from the Paper
#Figure 1A
#Figure 1C
#Figure 1D
#Figure 1F
#Figure 2A
#Figure 5A
#SFigure 8A

#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(BiocManager)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(scales)
library(data.table)
library(pathfindR)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library(vegan)
library(themetagenomics)
library(tidyr)
library(tableone)
library(stringr)

#Set Theme
theme<-
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), 
    legend.position="none")

#Choose Alpha/FDR
alpha = 0.01

##Load the files needed
file = "Merged.otu_table.biom"
map = "Map.COPD.SmNV.b2.txt"

# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))

#Merge abundance and meta data into a phyloseq object
lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)

#Make the columnnames of the pyloseq object of the phylogenetic tree
colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Load the tree file (use the unannotated.tree)
treefile = "97_otus.tree"
tree.obj = import_qiime(treefilename = treefile) 

#Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)

# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)

# # Create Genus Table
Genus.table = tax_glom(otu.table, taxrank = "Genus")

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}

#Create OTU and Genus Relative Tables
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)
Genus.rel.table = transformSampleCounts(Genus.table, normalizeSample)


#Subset Samples Genus
MetaTranscript.Genus.rel.table = subset_samples(Genus.rel.table, Subset_1==1)
MetaTranscript.Genus.rel.table = subset_samples(MetaTranscript.Genus.rel.table, Subject_Type_code %in% c(1, 2))
MetaTranscript.Genus.rel.table = subset_samples(MetaTranscript.Genus.rel.table, Metatranscriptome_plus_BKG  %in% c(1))
#Subset Samples OTU
MetaTranscript.otu.rel.table = subset_samples(otu.relative.table, Subset_1==1)
MetaTranscript.otu.rel.table = subset_samples(MetaTranscript.otu.rel.table, Subject_Type_code %in% c(1, 2))
MetaTranscript.otu.rel.table = subset_samples(MetaTranscript.otu.rel.table, Metatranscriptome_plus_BKG  %in% c(1))


#select most abundant taxa genera present in >2% relative abundance in 1% of the samples (this approach brings in > 70% of the data in almost all samples)
Genus.rel.wh1 = genefilter_sample(MetaTranscript.Genus.rel.table, filterfun_sample(function(x) x > 0.03), A = 0.01 * nsamples(MetaTranscript.Genus.rel.table))
Genus.rel.table1B = prune_taxa(Genus.rel.wh1, MetaTranscript.Genus.rel.table)
#set data tables  
GenusData <-otu_table(Genus.rel.table1B) #pruned to selected Genuses based on abundance

#cluster Genuses(row)
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")

#cluster samples(Col)
Samples.Bray.dist = distance(GenusData, method="bray")


#Set Color Scale for Heatmap
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))(100)
#Set Colors for each sample type for HeatMap
annon_colors= list(Sample_Type_Simple=c(BKG="#000001", BAL="#8EFA00", Sup="#932CE7"))

#Choose lables for Samples
df2 <- data.frame(Sample_Type_Simple = sample_data(Genus.rel.table1B)[,c("Sample_Type_Simple")], row.names = rownames(sample_data(Genus.rel.table1B)))

#Create dataframe of count data
df <- as.data.frame(GenusData)
#Get Taxa Names from Phyloseq Object
df = cbind(as(df, "data.frame"), as(tax_table(Genus.rel.table1B)[rownames(df), ], "matrix"))

#Replace any no genus annotation as NA
df[df=="g__"]<-NA
df[df=="f__"]<-NA
df[df=="o__"]<-NA
df[df=="c__"]<-NA
#Create name with family and (u.g)
df$gs <- ifelse(is.na(df$Genus),paste(df$Family,"(u.g.)"), paste(df$Genus))
df$gs <- ifelse(is.na(df$Family), paste(df$Order,"(u.g.)"),df$gs)
df$gs <- ifelse(is.na(df$Order), paste(df$Class,"(u.g.)"),df$gs)
df$gs <- ifelse(is.na(df$Class), paste(df$Phylum,"(u.g.)"),df$gs)

#Set Rownames
rownames(df) <- df$gs

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species")
df <- df[ , !(names(df) %in% drops)]
#Change the names if you need to
colnames(df) <- gsub("COPD.","",colnames(df))
colnames(df) <- gsub("SmNV.","",colnames(df))
colnames(df) <- gsub(".171","",colnames(df))
colnames(df) <- gsub(".172","",colnames(df))
colnames(df) <- gsub(".173","",colnames(df))

rownames(df2) <- gsub("COPD.","",rownames(df2))
rownames(df2) <- gsub("SmNV.","",rownames(df2))
rownames(df2) <- gsub(".171","",rownames(df2))
rownames(df2) <- gsub(".172","",rownames(df2))
rownames(df2) <- gsub(".173","",rownames(df2))

#---------
#Figure 1A
#---------
pdf("16S_Heatmap_Pruned.pdf", height = 15, width = 15)
    pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE,annotation_col=df2,scale="row",
    clustering_distance_rows = GenusData.Bray.dist,clustering_distance_cols = Samples.Bray.dist,
    clustering_method="average",
    gaps_col=50,
    border_color="black",
    color = colorRampPalette(c('#4169E1','#ffffff','#0000CD'))(100),
    annotation_colors=annon_colors[1],legend=FALSE)
dev.off()


#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(otu.table, ~ Sample_Type_DMM_Class_ReSeq)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#diagdds = estimateDispersions(diagdds)
#diagvst = getVarianceStabilizedData(diagdds)


#Subset BAL for analysis
diagdds <- diagdds[, diagdds$Subset_1 %in% c(1)]
diagdds <- diagdds[, diagdds$Subject_Type_code  %in% c(1, 2)]
diagdds <- diagdds[, diagdds$Metatranscriptome_plus_BKG  %in% c(1)]
diagdds.bal <- diagdds[, diagdds$Sample_Description_s_code  %in% c(5)]


write.table(assay(diagdds),file=  "16S_count_Data.txt", sep="\t", col.names = NA, row.names = TRUE)
write.table(colData(diagdds),file=  "16S_Mapping_File.txt", sep="\t", col.names = NA, row.names = TRUE)


#Make sure all unwanted levels are removed from dataset
diagdds$Sample_Type_DMM_Class_ReSeq <- droplevels(diagdds$Sample_Type_DMM_Class_ReSeq)


#Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
diagdds<- DESeq(diagdds)
res4 <- results(diagdds)
res4 <- results(diagdds, independentFiltering=FALSE)
#res4 <- results(dds4,cooksCutoff=TRUE)
#res4 <- results(dds4,pAdjustMethod="bonferroni")


# Reorder Results based on FDR
res4 = res4[order(res4$padj, na.last = NA), ]

#Create list of top 50 Significant Genes
select_genes = rownames(res4[res4$padj < alpha & !is.na(res4$padj), ])[1:50]

#Normalize Expression Data
exp  <- rlog(diagdds, fitType="local")
vsd <- varianceStabilizingTransformation(diagdds, blind=TRUE, fitType="local")


#=========================================================
////////////////////Alpha Diversity///////////////////////
#=========================================================#Calcultes Shannon Diversity
diagdds$Shannon = diversity(assay(diagdds), index = "shannon", MARGIN = 2, base = exp(1))
#Convert to data frame for ggplot
shannon = as.data.frame(colData(diagdds))
#Remove any zero values
shannon[shannon==0] <- NA

#Set Order Of Figure
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BKG", 1,NA)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BAL.BPT",2 ,shannon$or)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="BAL.SPT",3 ,shannon$or)
shannon$or <-ifelse(shannon$Sample_Type_DMM_Class_ReSeq=="Sup",4 ,shannon$or)

#---------
#Figure 1B
#---------
pdf("16S_Sample_Type_DMM_Class_SHANNON.pdf", height = 7, width = 5)
    ggplot(shannon, aes(x= reorder(Sample_Type_DMM_Class_ReSeq, +or), y=Shannon, color=Sample_Type_DMM_Class_ReSeq)) + 
    stat_boxplot(geom ='errorbar', width=0.1)+
    geom_boxplot(aes(color=Sample_Type_DMM_Class_ReSeq),outlier.shape = NA, width=0.5)+
    #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
    geom_jitter(aes(fill=Sample_Type_DMM_Class_ReSeq), position=position_jitter(0.2))+
    scale_fill_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) +
    scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) +
    scale_x_discrete(labels = c('BKG','BPT','SPT','UA'))+ 
    #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("Shannon Diversity") + 
    xlab("Sample Type")+
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.title=element_text(size=20,face="bold"),
    axis.text.x=element_text(colour = "black", size = rel(2)),
    axis.text.y=element_text(colour = "black", size = rel(2)),
    axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
    legend.position="none")
dev.off()

#Check Statistics
kruskal.test(Shannon ~ Sample_Type_DMM_Class_ReSeq, data = shannon)
Kruskal-Wallis chi-squared = 28.27, df = 3, p-value = 3.187e-06

#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================
#Create Distance Matrix
vegdist   = vegdist(t(assay(diagdds)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(diagdds), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Sample_Type_DMM_Class_ReSeq,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_DMM_Class_ReSeq",suffixes=c("",".centroid"))


#---------
#Figure 1D
#---------
pdf("16S_Sample_Type_DMM_Class_BRAY.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Sample_Type_DMM_Class_ReSeq)) +
    geom_point(size=5,alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Sample_Type_DMM_Class_ReSeq), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Sample_Type_DMM_Class_ReSeq))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.171", "COPD.0030.BAL.L.171", "COPD.0035.BAL.L.171") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL.BPT", "BAL.SPT", "BKG", "UA")), size=10) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),
    axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
    plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Check Statistics
adonis(vegdist ~ diagdds$Sample_Type_DMM_Class_ReSeq)
                                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
diagdds$Sample_Type_DMM_Class_ReSeq  3    5.4043 1.80145  5.8078 0.25464  0.001
Residuals                           51   15.8192 0.31018         0.74536
Total                               54   21.2235                 1.00000

diagdds$Sample_Type_DMM_Class_ReSeq ***
Residuals
Total
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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


#Create the tables for the function
Sup.Sample.otu.relative.table = subset_samples(MetaTranscript.otu.rel.table, Sample_Type_DMM_Class_ReSeq !="BKG")
Bray.dist = distance(Sup.Sample.otu.relative.table, method='bray')

######## Here, you set your parameters. In this case we use:   Bray  ###################
{
#set the sample data 
sample.data = sample_data(Sup.Sample.otu.relative.table)
#set the variable you want to do groupings by
Variable.Intergroup = sample_data(Sup.Sample.otu.relative.table)$Sample_Type_DMM_Class_ReSeq
#set the distance.matrix
distance.matrix = Bray.dist
##Select if you would like to add any other variables to the final table for clarifiation and processing. If none, then just repeat the variable.intergroup
#in this care s we need to add the Subject_ID, in this case: ID_Sample_Type_Subject_Type_Simple 
extraVar.addtotable = sample_data(Sup.Sample.otu.relative.table)$Sample_Type_DMM_Class_ReSeq
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


BAL.Sample.otu.relative.table = subset_samples(Sup.Sample.otu.relative.table, Sample_Type_DMM_Class_ReSeq %in% c("BAL.BPT","BAL.SPT"))
data.sup <- data.frame(data,Sample_Type=sample_data(BAL.Sample.otu.relative.table)$Sample_Type_DMM_Class_ReSeq,
                    Acetate    =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Acetate)),
                    Propionate =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Propionate)),
                    Butyrate   =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Butyrate)),
                    Isovalerate=as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Isovalerate)),
                    Valerate   =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Valerate)),
                    Hexanoate  =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Hexanoate)),
                    Octanoate  =as.numeric(as.character(sample_data(BAL.Sample.otu.relative.table)$Octanoate))
                    )

#Write Table For Figure 1I
write.table(data.sup, file = "Sup.BAL.Bray.txt", sep = "\t", row.names = FALSE)


------------
#BAL Only
------------
#Subset OTU Table
MetaTranscript.BAL.otu.table = subset_samples(otu.table, Subset_1==1)
MetaTranscript.BAL.otu.table = subset_samples(MetaTranscript.BAL.otu.table, Subject_Type_code %in% c(1, 2))
MetaTranscript.BAL.otu.table = subset_samples(MetaTranscript.BAL.otu.table, Metatranscriptome_plus_BKG  %in% c(1))
MetaTranscript.BAL.otu.table = subset_samples(MetaTranscript.BAL.otu.table, Sample_Description_s_code  %in% c(5))
MetaTranscript.BAL.rel.otu.table = transformSampleCounts(MetaTranscript.BAL.otu.table, normalizeSample)

out <- rarecurve(t(otu_table(MetaTranscript.BAL.otu.table)), step=50, col = "blue", cex = 0.6)

rare <- lapply(out, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- rownames(t(otu_table(MetaTranscript.BAL.otu.table)))

rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

pdf("16S_BAL_Rarefaction_Curves.pdf", height = 10, width = 10)
ggplot(data = rare)+
  geom_line(aes(x = raw.read, y = OTU, color = sample))+
  scale_x_continuous(labels =  scales::scientific_format())+
  xlab("Raw Reads")+
  ylab("OTU")+
  theme
dev.off()


#Subset BAL
diagdds.bal <- diagdds[, diagdds$Sample_Description_s_code  %in% c(5)]

#Make sure all unwanted levels are removed from dataset
diagdds.bal$Sample_Type_DMM_Class_ReSeq <- droplevels(diagdds.bal$Sample_Type_DMM_Class_ReSeq)

#Set Reference Level for Anlaysis
diagdds.bal$Sample_Type_DMM_Class_ReSeq <- relevel(diagdds.bal$Sample_Type_DMM_Class_ReSeq, ref ="BAL.BPT")

#Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
diagdds.bal<- DESeq(diagdds.bal, test="Wald", fitType="local")
res.bal <- results(diagdds.bal, cooksCutoff = FALSE)

# Reorder Results based on FDR
res.bal = res.bal[order(res.bal$padj, na.last = NA), ]

#Normalize Expression Data
exp.bal  <- rlog(diagdds.bal, fitType="local")
vsd.bal <- varianceStabilizingTransformation(diagdds.bal, blind=TRUE, fitType="local")


#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================
#Create Distance Matrix
vegdist   = vegdist(t(assay(diagdds.bal)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(diagdds.bal), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Sample_Type_DMM_Class_ReSeq,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type_DMM_Class_ReSeq",suffixes=c("",".centroid"))

#---------------------
#SFigure 8A
#---------------------
pdf("16S_BAL_Type_DMM_Class_BRAY_Propionate.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Sample_Type_DMM_Class_ReSeq)) +
    geom_point(aes(size=ifelse(is.na(as.numeric(as.character(Propionate))),1,as.numeric(as.character(Propionate))*5)),alpha=0.7) +
    scale_size_continuous(range=c(1,50))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#296218", "#EA3323")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Sample_Type_DMM_Class_ReSeq), size=0) +
    #geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Sample_Type_DMM_Class_ReSeq))+ 
    #labels centroids 
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL.BPT", "BAL.SPT")), size=10) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#-----------------------------------------------------
##Correlation Analysis with SCFA and Diversity Metrics
#-----------------------------------------------------

#Statistics with Propionate
adonis(vegdist ~ as.numeric(as.character(Propionate)), data.adonis)
                                     Df SumsOfSqs MeanSqs F.Model      R2
as.numeric(as.character(Propionate))  1    0.4983 0.49828  1.5827 0.08517
Residuals                            17    5.3521 0.31483         0.91483
Total                                18    5.8504                 1.00000
                                     Pr(>F)
as.numeric(as.character(Propionate))  0.067 .
Residuals
Total

#Spearman Beta
cor.test( ~ vegdist + as.numeric(as.character(Propionate)), 
         data=data.adonis ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

#Adonis Alpha
adonis(Shannon ~ as.numeric(as.character(Propionate)), shannon, na.rm=TRUE)
                                     Df SumsOfSqs  MeanSqs F.Model      R2
as.numeric(as.character(Propionate))  1   0.04662 0.046624 0.73434 0.04141
Residuals                            17   1.07934 0.063490         0.95859
Total                                18   1.12596                  1.00000
                                     Pr(>F)
as.numeric(as.character(Propionate))  0.098 .
Residuals
Total

#Spearman Alpha
cor.test( ~ Shannon + as.numeric(as.character(Propionate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
S = 1280, p-value = 0.6155
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
-0.122807

#Statistics with Acetate
#Adonis
data.adonis <- data.frame(colData(diagdds.bal))
adonis(vegdist ~ as.numeric(as.character(Acetate)), data.adonis)
                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Acetate))  1    0.5955 0.59553  1.9266 0.10179  0.055 .
Residuals                         17    5.2549 0.30911         0.89821 
Total                             18    5.8504                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Alpha
Shannon = diversity(assay(diagdds.bal), index = "shannon", MARGIN = 2, base = exp(1))
#Convert to data frame for ggplot
shannon = data.frame(colData(diagdds.bal))
#Adonis
adonis(Shannon ~ as.numeric(as.character(Acetate)), shannon, na.rm=TRUE)
                                  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Acetate))  1   0.06039 0.060393  0.9635 0.05364  0.075
Residuals                         17   1.06557 0.062680         0.94636
Total                             18   1.12596                  1.00000

as.numeric(as.character(Acetate)) .
Residuals

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Acetate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

S = 1242, p-value = 0.7155
alternative hypothesis: true rho is not equal to 0
sample estimates:
        rho
-0.08947368


#Statistics for Isovalerate
adonis(vegdist ~ as.numeric(as.character(Isovalerate)), data.adonis)
                                      Df SumsOfSqs MeanSqs F.Model      R2
as.numeric(as.character(Isovalerate))  1    0.5723 0.57233  1.8434 0.09783
Residuals                             17    5.2781 0.31047         0.90217
Total                                 18    5.8504                 1.00000
                                      Pr(>F)
as.numeric(as.character(Isovalerate))  0.068 .
Residuals
Total

#Adonis
adonis(Shannon ~ as.numeric(as.character(Isovalerate)), shannon, na.rm=TRUE)
                                      Df SumsOfSqs  MeanSqs F.Model      R2
as.numeric(as.character(Isovalerate))  1   0.05471 0.054712 0.86824 0.04859
Residuals                             17   1.07125 0.063015         0.95141
Total                                 18   1.12596                  1.00000
                                      Pr(>F)
as.numeric(as.character(Isovalerate))  0.112
Residuals
Total

#Spearman
cor.test( ~ Shannon + as.numeric(as.character(Isovalerate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

S = 1466, p-value = 0.2345
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.2859649

#Statistics with Butyrate
adonis(vegdist ~ as.numeric(as.character(Butyrate)), data.adonis)
                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Butyrate))  1    0.5378 0.53779  1.7209 0.09192  0.084
Residuals                          17    5.3126 0.31251         0.90808
Total                              18    5.8504                 1.00000

#Adonis
adonis(Shannon ~ as.numeric(as.character(Butyrate)), shannon, na.rm=TRUE)
                                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
as.numeric(as.character(Butyrate))  1    0.0878 0.087795  1.4377 0.07797  0.122
Residuals                          17    1.0382 0.061069         0.92203
Total                              18    1.1260                  1.00000

#SPEARMAN
cor.test( ~ Shannon + as.numeric(as.character(Butyrate)), 
         data=shannon ,
         method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

S = 1362, p-value = 0.4227
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
-0.1947368




#=========================================================
//////////////////////TABLES/////////////////////////////
#=========================================================
#Get Taxa Names from Phyloseq Object
res.bal = cbind(as(res.bal, "data.frame"), as(tax_table(MetaTranscript.BAL.otu.table)[rownames(res.bal), ], "matrix"))
#res.bal = cbind(as(res.bal, "data.frame"), as(tax_table(Species.bal.table)[rownames(res.bal), ], "matrix"))

#Replace OTU with Taxa
res.bal$row2 <- paste(res.bal$Domain,res.bal$Phylum,res.bal$Class,res.bal$Order,res.bal$Family,res.bal$Genus)
res.bal$row2 <- paste(res.bal$Domain,res.bal$Phylum,res.bal$Class,res.bal$Order,res.bal$Family,res.bal$Genus,res.bal$OTU)
#Replace Spaces with .
res.bal$row2 <- gsub('\\s+', '.', res.bal$row2)


#convert counts to a data.frame
cyber <- as.data.frame(assay(exp))

#Calculate Average Gene Expression per Condition
avg <- sapply(levels(dds5$SPT_BPT), function(lvl) rowMeans( counts(dds5,normalized=TRUE)[,dds5$SPT_BPT == lvl] ) )

#Convert Resuts table into a data.frame
res.bal <- as.data.frame(res.bal)

#Set Names of Results Table
res.bal <- setNames(cbind(rownames(res.bal), res.bal, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res.bal <- setNames(cbind(rownames(res.bal), res.bal, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","Species","row2"))

res.bal$names <- res.bal$Gene.symbol
res.bal$Gene.symbol <- res.bal$row2
res.bal$gs <- paste(res.bal$Genus,res.bal$OTU)
res.bal$gs <- gsub('\\s+', '.', res.bal$gs)

##add the OTU nunmber to the label
#Replace OTU with Taxa
res.bal$row3 <- paste(res.bal$Domain,res.bal$Phylum,res.bal$Class,res.bal$Order,res.bal$Family,res.bal$Genus, res.bal$OTU)
#Replace Spaces with .
res.bal$row3 <- gsub('\\s+', '.', res.bal$row3)

#rename
res.bal$Gene.symbol <- res.bal$row3




######get abundance data - mean relative - use otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE 
{
#decide what otu to save 
otu.to.save <-as.character(res.bal$names)

#from relative table we should get the mean across the row of the otu table
BAL.OTU.rel.table.df <- data.frame(otu_table(MetaTranscript.BAL.rel.otu.table))
BAL.OTU.rel.table.df.meanRA <- rowMeans(BAL.OTU.rel.table.df)

#need to subset AND reorder just the otus that we have 
BAL.OTU.rel.table.df.meanRA.save <- BAL.OTU.rel.table.df.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res.bal$abundance <- BAL.OTU.rel.table.df.meanRA.save
}



keepOTUs = rownames(res.bal[res.bal$padj < alpha & !is.na(res.bal$padj), ])
keepOTUs = res.bal[res.bal$adj.P.Val < alpha & res.bal$logFC>2, "gs"]
keepOTUs = res.bal[res.bal$adj.P.Val < alpha & res.bal$logFC>2, "Genus"]

#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================
# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res.bal$sig <- -log10(res.bal$adj.P.Val)
sum(is.infinite(res.bal$sig))

res.bal[is.infinite(res.bal$sig),"sig"] <- 350

# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(res.bal$pvalue)

# sum(genes.to.plot)
range(res.bal[genes.to.plot, "logFC"])

## Volcano plot of adjusted p-values
cols <- densCols(res.bal$logFC, res.bal$sig)
cols[res.bal$pvalue ==0] <- "purple"
cols[res.bal$logFC > 2 & res.bal$adj.P.Val < alpha ] <- "red"
res.bal$pch <- 19
res.bal$pch[res.bal$pvalue ==0] <- 6

gn.selected <- abs(res.bal$logFC) > 10 & res.bal$adj.P.Val < alpha 

#----------
#Figure 2A
#----------
pdf(file="16S_BAL_SPT_vs_BPT_Volcano_Plot_rel_Abundance_Small_FDR_0.05.pdf", width=5, height=5)
    ggplot(res.bal, aes(x = logFC, y = sig,label=Gene.symbol)) +
    geom_point(color=cols, size = ifelse(res.bal$adj.P.Val < alpha, 200 * res.bal$abundance, 2), alpha=0.5) + #Chose Colors and size for dots
    geom_text_repel(aes(label=ifelse(res.bal$logFC>8 & res.bal$adj.P.Val < alpha & res.bal$Genus!="g__", as.character(res.bal$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
    xlab("Effect size: log2(fold-change)") + #label X Axis
    ylab("-log10(adjusted p-value)") + #label Y Axis
    ylim(0,20)+
    theme #Set Theme
dev.off() 



#////////////////////////////////////////////////////
#------------------------------------>PICRUST Analysis
#////////////////////////////////////////////////////
#Not the Relative Abundance Table
otu.table <- t(as.matrix((otu_table(MetaTranscript.BAL.otu.table))))

#Download the Reference KO Table
tmp <- tempdir()
download_ref(tmp,reference='gg_ko',overwrite=FALSE)

#Run the PICRUST analysis with the OTU Table and the Reference KO Table
#system.time to see how long it takes
system.time(FUNCTIONS <- picrust(otu.table,rows_are_taxa=FALSE,
                                 reference='gg_ko',reference_path=tmp,
                                 cn_normalize=FALSE,sample_normalize=FALSE,
                                 drop=TRUE))
#KO Annotation
KO <- as.data.frame(t(FUNCTIONS$fxn_table))
#KEGG Pathway
pathway <- FUNCTIONS$fxn_meta$KEGG_Pathways
#Add pathway to table
KO$path <- sapply(pathway, paste0, collapse=",") 
#Remove Special Characters
KO$path <- gsub("),c\\(","|",KO$path)
KO$path <- gsub(",",";",KO$path)
KO$path <- gsub("c\\(","",KO$path)
KO$path <- gsub(")","",KO$path)
KO$path <- gsub("[\"]","",KO$path)
#Remove White Space
library(stringr)
KO$path  <- str_squish(KO$path)
#add KEGG ID to the pathway
KO$path  <- paste0(rownames(KO),"|",KO$path)
#change Rownames to full pathway
rownames(KO) <- KO$path
#remove path column
KO <- KO[,!names(KO) %in% "path"]

#Load Meta Data
coldata <- read.delim2("Map.COPD.SmNV.b2.txt", sep="\t")

#Subset only samples of interest
coldata = coldata[coldata$Subset_1==1,]
coldata = coldata[coldata$Subject_Type_code %in% c(1, 2),]
coldata = coldata[coldata$Metatranscriptome_plus_BKG  %in% c(1),]
#coldata = coldata[coldata$Sample_Description_s_code  %in% c(5),]

#Order Meta Data by SampleId
coldata <- coldata[order(coldata$X),]

#Order Count Data by SampleID
KO <-KO[, order(colnames(KO))]

#Confirm Sample IDs match for Count and Meta Data
table(colnames(KO)==as.character(coldata$X))

#Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds <- DESeqDataSetFromMatrix(countData = KO,
                              colData = coldata,
                              design= ~ Sample_Type_DMM_Class_ReSeq)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds), 1, gm_mean)

#Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset
dds$Sample_Type_DMM_Class_ReSeq <- droplevels(dds$Sample_Type_DMM_Class_ReSeq)

#Set Reference Level for Anlaysis
dds$Sample_Type_DMM_Class_ReSeq <- relevel(dds$Sample_Type_DMM_Class_ReSeq, ref ="BAL.BPT")

#Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
dds<- DESeq(dds)
res <- results(dds, cooksCutoff = FALSE)

#Run Transformation
rld <- rlog(dds, fitType="local")
rld30 <- ifelse(assay(rld)<0,0,assay(rld))

#Select top KOs
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

#Convert Results into a dataframe
results <- as.data.frame(res)
#Create Column of Rownames
results$KO_Number <- rownames(results)
#Keep only the KO Number
results$KO_Number <- substr(results$KO_Number, 0, 6)
#Asign Hierarchy to KO Number
results <- assign_hierarchy(count_data=results, keep_unknowns=TRUE, identifier="KO_Number")
#Remove Number from Subclass
results$KO_Subclass_2 <- gsub('[[:digit:]]+', '', results$KO_Subclass_2 )
#Remove Path Number from Subclass
results$KO_Subclass_2 <- gsub("\\s*\\[[^\\)]+\\]","",results$KO_Subclass_2)


#=========================================================
//////////////////////TABLES/////////////////////////////
#=========================================================
#Convert Resuts table into a data.frame
res <- as.data.frame(results)

######get abundance data - mean relative - use otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE 
{
#decide what otu to save 
ko.to.save <-as.character(rownames(res))

#from relative table we should get the mean across the row of the otu table
BAL.ko.table.df <- data.frame(assay(dds))
BAL.ko.table.df.meanRA <- rowMeans(BAL.ko.table.df)

#need to subset AND reorder just the otus that we have 
BAL.ko.table.df.meanRA.save <- BAL.ko.table.df.meanRA[ko.to.save]

#add the abundnace data for the res dataframe
res$abundance <- BAL.ko.table.df.meanRA.save
}

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "KO_Number","KO_Class", "KO_Subclass_1","KO_Subclass_2","abundance")) 
res$path <- paste(res$KO_Number,res$KO_Class,res$KO_Subclass_1,res$KO_Subclass_2, sep="_")

#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================
# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res$sig <- -log10(res$adj.P.Val)
sum(is.infinite(res$sig))

res[is.infinite(res$sig),"sig"] <- 350

## Volcano plot of adjusted p-values
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

#----------
#Figure 1F
#----------
pdf(file="PICRUST_BAL_SPT_vs_BPT_Volcano_Plot_Small_FDR_0.05.pdf", width=5, height=5)
    ggplot(res, aes(x = logFC, y = sig,label=path)) +
    geom_point(color=cols, size = ifelse(res$adj.P.Val < alpha, 10^-3.7*res$abundance, 2), alpha=0.5) + #Chose Colors and size for dots
    geom_text_repel(aes(label=ifelse(res$logFC>0 & res$adj.P.Val < 0.000000000000000001 & !is.na(res$KO_Subclass_2) , as.character(res$path),'')),size=2,force=25, segment.colour="grey",segment.alpha=0.5) +
    theme(legend.position = "none") +
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
    xlab("Effect size: log2(fold-change)") +
    ylab("-log10(adjusted p-value)") + 
    ylim(0,20)+
    theme
dev.off() 



#------------------------------------------------
#---->Using SCFA analysis with Inferred MetaGenome data
#------------------------------------------------
#Get KEGG Table from Picrust
kegg <- KO
#Create Column for KO
kegg$ko <- substr(rownames(kegg), 1,6)

#Fix Rownames of database
rownames(kegg) <- kegg$ko
#Drop Reddundant Columns
drops <- c("ko")
kegg <- kegg[ , !(names(kegg) %in% drops)]

#Make CountData and MetaData into DESEq Object; Choose the comparison Variable as design
dds.keggscfa <- DESeqDataSetFromMatrix(countData = kegg,
                              colData = coldata,
                              design= ~ Sample_Type_DMM_Class_ReSeq)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds.keggscfa), 1, gm_mean)

#Estimate Factors of DESeq Object
dds.keggscfa <- estimateSizeFactors(dds.keggscfa, geoMeans = geoMeans)

#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds.keggscfa, normalized=TRUE) >= 100 ) >= 3
dds.keggscfa <- dds.keggscfa[idx,]

#=========================================================
//////////////DOT PLOTS OF SIGNFICANT GENES////////////////
#=========================================================
#Convert Normalized Assay Data into a Data Table
as <- data.table(melt(assay(dds.keggscfa)))

#Create a Variable name to match with MetaData
as$SampleID <- as$Var2

#Create a Data Frame of MetaData
assay <- as.data.frame(colData(dds.keggscfa))
assay$SampleID <- rownames(assay)

#Merge Assay Data and Meta Data
as2 <- merge(as,assay, all=TRUE)

#convert Genus to a character vector from a factor
as2$KEGG <- as.character(as2$Var1)

#Create Data Frame copy
as3 <- as2

#Extract KEGGs for Acetate, Propionate and Butyrate
keepKOSC <- c("K01738","K00925","K01034")
as2 <- as3[as3$KEGG %in% keepKOSC,]

#Test Overall Significance
kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2) 
Kruskal-Wallis chi-squared = 10.097, df = 3, p-value = 0.01776

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2.BPT.SPT)
W = 285.5, p-value = 0.05718

#Replace 0 with 0.1
as2$value <- ifelse(as2$value==0,0.1,as2$value)
#Change BKG label so that it comes first in order
as2$Sample_Type_DMM_Class_ReSeq <- as.character(as2$Sample_Type_DMM_Class_ReSeq)
as2$Sample_Type_DMM_Class_ReSeq <- ifelse(as2$Sample_Type_DMM_Class_ReSeq=="BKG","aBKG",as2$Sample_Type_DMM_Class_ReSeq)
#Change Order so its Acetate Propionate and Butyrate
as2$order <-ifelse(as2$KEGG=="K01738",1,
            ifelse(as2$KEGG=="K00925",2,3))

#---------
#Figure 5B
#---------
pdf("PICRUST_KEGG_available_BAL_SPT_vs_BPT_BOX_Plot.pdf.pdf", height = 7, width = 10)
    ggplot(as2, aes(x=reorder(KEGG,+order), y=value, color=Sample_Type_DMM_Class_ReSeq)) +
    geom_boxplot(aes(color=Sample_Type_DMM_Class_ReSeq),outlier.shape = NA, width=0.7)+
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

kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2) 
Kruskal-Wallis chi-squared = 16.653, df = 3, p-value = 0.0008329

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2.BPT.SPT)
W = 17, p-value = 0.02202

#Extract Propionate Data
keepKOSC <- c("K00925")
as2 <- as3[as3$KEGG %in% keepKOSC,]

kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2) 
Kruskal-Wallis chi-squared = 20.909, df = 3, p-value = 0.0001099

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2.BPT.SPT)
W = 10, p-value = 0.002988

#Extract Butyrate Data
keepKOSC <- c("K01034")
as2 <- as3[as3$KEGG %in% keepKOSC,]

kruskal.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2) 
Kruskal-Wallis chi-squared = 8.2221, df = 3, p-value = 0.04164

as2.BPT.SPT <- as2[as2$Sample_Type_DMM_Class_ReSeq %in% c("BAL.BPT", "BAL.SPT"),]

wilcox.test(value ~ Sample_Type_DMM_Class_ReSeq, data = as2.BPT.SPT)
W = 23.5, p-value = 0.08627



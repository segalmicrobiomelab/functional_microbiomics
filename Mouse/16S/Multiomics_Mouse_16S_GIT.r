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
library(vegan)
library(dplyr)
library(MetaboSignal)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(htmlwidgets)
library(plotly)
library(picante)
library(phyloseq)
library(cowplot)
library(grid)
library(gridExtra) 

#Set Theme for Figures
theme<-	theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.05

#/////////////////////////////////////////////////////
#---------------------------------->16S Analysis
#/////////////////////////////////////////////////////
#Set Working directory
setwd("/gpfs/data/segallab/home/sulaii01/MultiOmics/Mouse/16S/")

#Load File
load(file="Mouse.16S.RData")

##Load the files needed
file = "otu_table.biom"
map = "MSQ.109.Map_corrected.txt"

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

#Only Mouse Experiment
otu.table = subset_samples(otu.table, Experiment %in% "Mouse")

#Make Genus Table
Genus.table = tax_glom(otu.table, taxrank = "Genus")


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}

#Relative Abundance Tables
otu.rel.table = transformSampleCounts(otu.table, normalizeSample)
Genus.rel.table = transformSampleCounts(Genus.table, normalizeSample)


#BAL Tables
otu.BAL.table = subset_samples(otu.table, Description %in% c("BAL","MOC"))
Genus.BAL.table = subset_samples(Genus.table, Description %in% c("BAL","MOC"))
Genus.BAL.PBSMOC.table = subset_samples(Genus.table, Time_Exp %in% c("1_Hour_PBS","Innoc_MOC","Innoc_PBS"))
Genus.BAL.MOC.table = subset_samples(Genus.table, Time_Exp %in% c("Innoc_MOC"))
otu.BAL.MOC.table = subset_samples(otu.table, Time_Exp %in% c("Innoc_MOC"))
otu.BAL.PBS.table = subset_samples(otu.table, Time_Exp %in% c("Innoc_PBS"))

otu.BAL.rel.table = subset_samples(otu.rel.table, Description %in% c("BAL","MOC"))

Genus.BAL.rel.table = subset_samples(Genus.rel.table, Description %in% c("BAL","MOC"))


#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds.bal <- phyloseq_to_deseq2(Genus.BAL.table, ~ Time_Exp)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds.bal), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds.bal = estimateSizeFactors(diagdds.bal, geoMeans = geoMeans)

#Subset BAL for analysis
diagdds <- diagdds[, diagdds$Subset_1 %in% c(1)]
diagdds <- diagdds[, diagdds$Subject_Type_code  %in% c(1, 2)]
diagdds <- diagdds[, diagdds$Metatranscriptome_plus_BKG  %in% c(1)]
diagdds.bal <- diagdds[, diagdds$Sample_Description_s_code  %in% c(5)]

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
centroids <- aggregate(cbind(PC1,PC2)~Time_Exp,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Time_Exp",suffixes=c("",".centroid"))

#--------------------
#-----------Figure 6B
#--------------------
pdf("16S_BAL_DMM_Class_BRAY.pdf", height = 7, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Time_Exp)) +
    geom_point(size=5,alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    #scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Time_Exp), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Time_Exp))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.171", "COPD.0030.BAL.L.171", "COPD.0035.BAL.L.171") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("1 Day", "1 Hour", "PBS", "3 Days","4 Hours","7 Days","MOC")), size=10) +
    scale_color_manual(values=c("#6ABD23","#C49A02","#FB83DF","#F8766D","#18C59D","#06B9EB","#A88EFF")) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_line(color="grey",size=0.2,linetype=3),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),
    axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
    plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

adonis(vegdist ~ diagdds.bal$Time_Exp)
                     Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
diagdds.bal$Time_Exp  7    4.1995 0.59992   7.324 0.7855  0.001 ***
Residuals            14    1.1468 0.08191         0.2145
Total                21    5.3462                 1.0000
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#=========================================================
///////////////BARPLOT OF TAXA////////////////
#=========================================================
# function to find the most abundant taxa
# Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
# and identifies which taxa is most abundant for which sample
#Find the top 3 most differentially enriched Taxa 
keepOTUs = names(sort(taxa_sums(otu.BAL.MOC.table), TRUE)[1:3])

#select most abundant taxa genera present in >2% relative abundance in 1% of the samples (this approach brings in > 70% of the data in almost all samples)
genus.wh1 = genefilter_sample(Genus.BAL.rel.table, filterfun_sample(function(x) x > 0.03), A = 0.01 * nsamples(Genus.BAL.rel.table))
genus.table1B = prune_taxa(genus.wh1, Genus.BAL.rel.table)

#OTU TABLe pruning 
otu.wh1 = genefilter_sample(otu.BAL.rel.table, filterfun_sample(function(x) x > 0.03), A = 0.01 * nsamples(otu.BAL.rel.table))
otu.table1B = prune_taxa(otu.wh1, otu.BAL.rel.table)

Genus.MOC.table1B <- subset_samples(genus.table1B, Time_Exp %in% c("Innoc_MOC"))
otu.MOC.table1B <- subset_samples(otu.table1B, Time_Exp %in% c("Innoc_MOC"))

keepOTUs = names(sort(taxa_sums(Genus.MOC.table1B), TRUE)[1:3])
keepOTUs = names(sort(taxa_sums(otu.MOC.table1B), TRUE)[1:3])



#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(genus.table1B))), ntaxa(genus.table1B)), genus.table1B)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

    x10 = prune_taxa(tail(names(sort(taxa_sums(otu.table1B))), ntaxa(otu.table1B)), otu.table1B)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

    x10 = prune_taxa(tail(names(sort(taxa_sums(Genus.LT.table))), ntaxa(Genus.LT.table)), Genus.LT.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comaring Status to a character vector from a factor
dat$Location <- as.character(dat$Time_Exp)


# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)

dat$Genus <- dat$o
dat <- dat[!is.na(dat$Genus),]


#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]
data <- setDT(dat)[,list(Abundance=as.numeric(mean(Abundance, na.rm=TRUE)), iqr=as.numeric(sd(Abundance, na.rm=TRUE))), by=c("Location", "Genus")]


#Remove MOC & PBS
data2 <- data[data$Location!=c("Innoc_MOC","1_Hour_PBS" ]

#Calculate the Other Taxa rel Abundance
data3 <-
    data2 %>%
    bind_rows(data2 %>% 
            group_by(Location) %>%
            summarise_if(is.numeric, funs(sum))) %>% 
    mutate(Genus = ifelse(is.na(Genus), "Other", Genus))

#Covert Other to 1- and covert all to Percentage
data3$Abundance <- ifelse(data3$Genus=="Other",(1-data3$Abundance)*100,data3$Abundance*100)
data3$or <-ifelse(data3$Location=="1_Hour_MOC",0 ,NA)
data3$or <-ifelse(data3$Location=="4_Hours_MOC",4 ,data3$or)
data3$or <-ifelse(data3$Location=="1_Day_MOC",12 ,data3$or)
data3$or <-ifelse(data3$Location=="3_Days_MOC",36 ,data3$or)
data3$or <-ifelse(data3$Location=="7_Days_MOC",84 ,data3$or)


#---------------
#-----Figure 6D
#---------------
pdf("16S_Relative_Abundance_Median_Stacked_BAL_OTU_over_Time.pdf", height = 7, width = 10)
    ggplot(data3, aes(x=or, y=Abundance)) +
    #geom_blank() +
    geom_bar(stat="identity",width=2.5,aes(fill=factor(Genus,
    levels=c("Other","f__Prevotellaceae_g__Prevotella_4307652","f__Streptococcaceae_g__Streptococcus_1083194",
    "f__Veillonellaceae_g__Veillonella_585419"))))+
    #facet_grid(. ~ Genus, scales="free",labeller=labeller(Genus=labels)) +
    #facet_wrap(. ~ Genus, scales="fixed",labeller=labeller(Genus=labels),ncol=1) +
    #geom_errorbar(aes(ymin=value, ymax=iqr,color=Genus), width=.2,position=position_dodge(.9)) +
    xlab("") +
    ylab("% Relative Abundance")+
    scale_y_continuous(limits=c(0, 105), expand = c(0, 0))+
    #scale_y_continuous(expand = c(0, 0))+
    #scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    scale_x_discrete(labels=c("Innoc_MOC"="MOC", "1_Hour_MOC"="1Hr", "4_Hours_MOC"= "4Hr",
                                "1_Day_MOC"="1D", "3_Days_MOC"="3D","7_Days_MOC"="7D","1_Hour_PBS"="PBS"))+
    scale_fill_manual(values=c( "#D6D6D6","#D6504C","#47B665","#FBBE51")) +
    #scale_color_manual(values=c("#D6504C","#47B665","#FBBE51")) +
    #coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size=20,face="bold"),panel.spacing.x=unit(2.5, "lines"),axis.title=element_text(size=20, face="bold"))
dev.off()

#MOC PIE CHART
MOC <- data[data$Location=="Innoc_MOC"]
#Calculate Other
MOC <-
    MOC %>%
    bind_rows(MOC %>% 
            summarise_if(is.numeric, funs(sum))) %>% 
    mutate(Genus = ifelse(is.na(Genus), "Other", Genus))

MOC$Abundance <- ifelse(MOC$Genus=="Other",(1-MOC$Abundance),MOC$Abundance)

#Calculate YMax and YMin for Pie Chart
MOC$ymax = cumsum(MOC$Abundance)
MOC$ymin = c(0, head(MOC$ymax, n=-1))

#---------------
#------Figure 6A
#---------------
pdf("16S_MOC_Donut.pdf", height = 5, width = 10)
    ggplot(MOC, aes(fill=Genus, ymax=ymax, ymin=ymin, xmax=8, xmin=3)) +
    geom_rect() +
    coord_polar(theta="y") +
    xlim(c(0, 8)) +
    scale_fill_manual(values=c("#D6504C","#47B665","#FBBE51","#D6D6D6"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.ticks= element_blank(),axis.text=element_blank(),legend.position="none")
dev.off()


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


library(ggpmisc)
library(dplyr)


#Create the tables for the function
Sup.Sample.otu.relative.table = subset_samples(MetaTranscript.otu.rel.table, Sample_Type_DMM_Class_ReSeq !="BKG")
Bray.dist = distance(Genus.BAL.rel.table, method='bray')

######## Here, you set your parameters. In this case we use:   Bray  ###################
{
#set the sample data 
sample.data = sample_data(Genus.BAL.rel.table)
#set the variable you want to do groupings by
Variable.Intergroup = sample_data(Genus.BAL.rel.table)$Time_Exp
#set the distance.matrix
distance.matrix = Bray.dist
##Select if you would like to add any other variables to the final table for clarifiation and processing. If none, then just repeat the variable.intergroup
#in this care s we need to add the Subject_ID, in this case: ID_Sample_Type_Subject_Type_Simple 
extraVar.addtotable = sample_data(Genus.BAL.rel.table)$Time_Exp
#set the file name you'd like
filename = "Unique.ID.Bronch.Cohort.Bray.txt"
allEncompassing.df = NULL


#run it for all distances between samples
intergroup.distances(sample.data, Variable.Intergroup, distance.matrix, extraVar.addtotable, filename)

#Remove any BAL to BAL comparissons
allEncompassing.df <- allEncompassing.df[allEncompassing.df$SampleType2 %in% c("1_Hour_PBS","Innoc_MOC"),]

#Create a Variable that is BAL to BKG Comparisson
allEncompassing.df$comparison <- ifelse(allEncompassing.df$Subject2=="Innoc_MOC",
                                paste(allEncompassing.df$Subject1), paste(allEncompassing.df$Subject2))
#Remove Duplicates
allEncompassing.df = allEncompassing.df[allEncompassing.df$comparison=="1_Hour_PBS",]

#Create a variable for SampleID
allEncompassing.df$SampleID <- ifelse(allEncompassing.df$Subject2=="1_Hour_PBS",
                                as.character(allEncompassing.df$Subject1), as.character(allEncompassing.df$Subject2))

#Keep only the Sample ID and the distance
allEncompassing.df <- allEncompassing.df[,names(allEncompassing.df) %in% c("SampleID","Distance")]

#Summarize the Mean of the data per sampleID
data <- allEncompassing.df %>%
  group_by(SampleID) %>%
  summarize(mean_size = mean(Distance, na.rm = TRUE), 
            sd= sd(Distance, na.rm=TRUE),
            median=median(Distance,na.rm=TRUE),
            iqr=quantile(Distance, probs=.75, na.rm=TRUE))

#Set Order Of Figure
data$or <-ifelse(data$SampleID=="Innoc_MOC",0 ,NA)
data$or <-ifelse(data$SampleID=="1_Hour_MOC",0 ,data$or)
data$or <-ifelse(data$SampleID=="4_Hours_MOC",2 ,data$or)
data$or <-ifelse(data$SampleID=="1_Day_MOC",12 ,data$or)
data$or <-ifelse(data$SampleID=="3_Days_MOC",36 ,data$or)
data$or <-ifelse(data$SampleID=="7_Days_MOC",84 ,data$or)

data$group <- ifelse(data$SampleID=="Innoc_MOC", "1.MOC","2.Others")

moc  <- data[data$SampleID=="Innoc_MOC",]
time <- data[data$SampleID!="Innoc_MOC",]

#Plot Median + IQR
    #ggplot(data, aes(x=reorder(SampleID,+or), y=mean_size, fill=SampleID)) +
mocplot <-    
    ggplot(moc, aes(x=or, y=mean_size, fill=SampleID)) +
    #geom_blank() +
    geom_bar(stat="identity", width =0.2)+
    #facet_grid(~group,scales="free",space="free") + 
    #facet_grid(~group,scales="free",space="free") + 
    #facet_grid(. ~ Genus, scales="free",labeller=labeller(Genus=labels)) +
    geom_errorbar(aes(ymin=mean_size, ymax=sd,color=SampleID), width=.2, position=position_dodge(.9)) +
    xlab("") +
    ylab("Mean Distance from PBS")+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_discrete(labels=c("MOC"))+
    scale_fill_manual(values=c("#A88EFF")) +
    scale_color_manual(values=c("#A88EFF")) +
    #coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size=20,face="bold"),panel.spacing.x=unit(2.5, "lines"),axis.title=element_text(size=20, face="bold"))


timeplot <-    #ggplot(data, aes(x=reorder(SampleID,+or), y=mean_size, fill=SampleID)) +
    ggplot(time, aes(x=or, y=mean_size, fill=SampleID)) +
    #geom_blank() +
    geom_area(stat = "identity", color= "black",fill="lightgrey", linetype="dashed",alpha=0.2)+
    geom_bar(stat="identity", width =1)+
    #facet_grid(~group,scales="free",space="free") + 
    #facet_grid(~group,scales="free",space="free") + 
    #facet_grid(. ~ Genus, scales="free",labeller=labeller(Genus=labels)) +
    geom_errorbar(aes(ymin=mean_size, ymax=sd,color=SampleID), width=.2, position=position_dodge(.9)) +
    xlab("") +
    ylab("")+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_discrete(labels=c("1_Hour_MOC"="1 Hour", "4_Hours_MOC"= "4 Hours",
                                "1_Day_MOC"="1 Day", "3_Days_MOC"="3 Days","7_Days_MOC"="7 Days"))+
    scale_fill_manual(values=c("#6ABD23","#C49A02","#F8766D","#18C59D","#06B9EB")) +
    scale_color_manual(values=c("#6ABD23","#C49A02","#F8766D","#18C59D","#06B9EB")) +
    #coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line.y = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.line.x = element_line(colour = "black"),axis.ticks.x=element_line(colour = "black"),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=unit(2.5, "lines"),
    axis.title=element_text(size=20, face="bold"))

#-----------------
#------Figure 6C
#-----------------
pdf("16S_InterGroup_Distance_From_PBS.pdf", height = 5, width = 10)
plot_grid(mocplot, timeplot, align = "h", ncol = 2, rel_widths = c(1.7/10, 8.3/10))
dev.off()



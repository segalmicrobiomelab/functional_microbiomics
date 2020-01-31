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
library(cowplot)
library(funrar)
library(ggpmisc)
library(dplyr)


#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.01

#Load Meta Data
coldata <- read.delim2("Map.Mouse2.txt", sep="\t")

#Remove Sample with 0 Reads for all Genes
coldata <- coldata[coldata$Mouse.number!="M18",]

#Remove Sample with dimer peak too large
coldata <- coldata[coldata$Mouse.number!="M15",]

#Order Meta Data by SampleId
coldata <- coldata[order(coldata$Mouse.number),]

#load Count Data
mycounts <-read.delim2("KEGG_gene_table.txt", sep="\t", row.names=1)

#Remove dimer peak too large
mycounts = mycounts[,!(names(mycounts) %in% "M15")]

#Order Count Data by SampleID
mycounts <-mycounts[, order(colnames(mycounts))]

#Confirm Sample IDs match for Count and Meta Data
table(colnames(mycounts)==as.character(coldata$Mouse.number))

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0

#Create Copy of Count Table
mycounts2 <- mycounts

#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts, function(x) as.numeric(as.character(x))),
                   check.names=F, row.names = rownames(mycounts3))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

dds <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata,
                              design= ~ Time_exp)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds), 1, gm_mean)

#Estimate Factors of DESeq Object
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)

#For genes with lower counts, however, the values are shrunken towards the genes' averages across all sample
rld <- rlog(dds, fitType="local")
vsd <- varianceStabilizingTransformation(dds)

#Make sure all unwanted levels are removed from dataset
dds$Time_exp <- droplevels(dds$Time_exp)
rld$Time_exp <- droplevels(rld$Time_exp)
vsd$Time_exp <- droplevels(vsd$Time_exp)

#keep Remved ungrouped KEGGS
dds1 <- dds[!grepl("UNGROUPED|UNMAPPED",rownames(assay(dds))),]
rld1 <- rld[!grepl("UNGROUPED|UNMAPPED",rownames(assay(rld))),]
vsd1 <- vsd[!grepl("UNGROUPED|UNMAPPED",rownames(assay(vsd))),]
#keep only data which includes Taxa Data
dds2 <- dds1[grep("g_|unclassified",rownames(assay(dds1))),]
rld2 <- rld1[grep("g_|unclassified",rownames(assay(rld1))),]
vsd2 <- vsd1[grep("g_|unclassified",rownames(assay(vsd1))),]
#keep only data for KEGG
dds3 <- dds1[!grepl("g_|unclassified",rownames(assay(dds1))),]
rld3 <- rld1[!grepl("g_|unclassified",rownames(assay(rld1))),]
vsd3 <- vsd1[!grepl("g_|unclassified",rownames(assay(vsd1))),]

#Remove Negative Controls for analysis
dds2noneg <- dds2[, dds2$Mouse.number!= "NEG2"]
dds2noneg <- dds2noneg[, dds2noneg$Mouse.number!= "NEG1"]

dds3noneg <- dds3[, dds3$Mouse.number!= "NEG2"]
dds3noneg <- dds3noneg[, dds3noneg$Mouse.number!= "NEG1"]

rld3noneg <- rld3[, rld3$Mouse.number!= "NEG2"]
rld3noneg <- rld3noneg[, rld3noneg$Mouse.number!= "NEG1"]


#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================
#Create Distance Matrix
vegdist   = vegdist(t(assay(rld3noneg)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(rld3noneg), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Time_exp,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Time_exp",suffixes=c("",".centroid"))

#---------------
#------Figure 6E
#---------------
pdf("KEGG_no_M15_Metatranscriptome_BAL_Time_BRAY.pdf", height = 7, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Time_exp)) +
    geom_point(size=5,alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    #scale_color_manual(values=c("#296218", "#EA3323", "#000000","#932CE7")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Time_exp), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Time_exp))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.171", "COPD.0030.BAL.L.171", "COPD.0035.BAL.L.171") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("1 Day", "1 Hour", "3 Days","4 Hours","7 Days","MOC", "PBS")), size=10) +
    scale_color_manual(values=c("#6ABD23","#C49A02","#F8766D","#18C59D","#06B9EB","#A88EFF","#FB83DF")) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Time_exp), size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_line(color="grey",size=0.2,linetype=3),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),
    axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
    plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
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
#######################   DON'T CHANGE THIS FUNCTION: STARTS HERE     ##################################################

#Create the tables for the function
mat <- ifelse(assay(rld3noneg)<0,0,assay(rld3noneg))
Bray.dist   = vegdist(t(mat), method="bray")

######## Here, you set your parameters. In this case we use:   Bray  ###################
{
#set the sample data 
sample.data = colData(rld3noneg)
#set the variable you want to do groupings by
Variable.Intergroup = colData(rld3noneg)$Time_exp
#set the distance.matrix
distance.matrix = Bray.dist
##Select if you would like to add any other variables to the final table for clarifiation and processing. If none, then just repeat the variable.intergroup
#in this care s we need to add the Subject_ID, in this case: ID_Sample_Type_Subject_Type_Simple 
extraVar.addtotable = colData(rld3noneg)$Time_exp
#set the file name you'd like
filename = "Unique.ID.Bronch.Cohort.Bray.txt"
allEncompassing.df = NULL


#run it for all distances between samples
intergroup.distances(sample.data, Variable.Intergroup, distance.matrix, extraVar.addtotable, filename)


#Create Comparisson Variable
allEncompassing.df$comparison <- ifelse(allEncompassing.df$Subject1=="PBS",
                                paste(allEncompassing.df$Subject2,"to",allEncompassing.df$Subject1), paste(allEncompassing.df$Subject1,"to",allEncompassing.df$Subject2))
#Keep only PBS Comparisson
allEncompassing.df <- allEncompassing.df[grepl("PBS",allEncompassing.df$comparison),]

#Create a variable for SampleID
allEncompassing.df$SampleID <- ifelse(allEncompassing.df$Subject2=="PBS",
                                as.character(allEncompassing.df$Subject1), as.character(allEncompassing.df$Subject2))

#Keep only the Sample ID and the distance
allEncompassing.df <- allEncompassing.df[,names(allEncompassing.df) %in% c("SampleID","Distance")]

#Summarize the Mean of the data per sampleID
data <- allEncompassing.df %>%
  group_by(SampleID) %>%
  summarize(mean_size = mean(Distance, na.rm = TRUE), sd= sd(Distance, na.rm=TRUE))

#Set Order Of Figure
data$or <-ifelse(data$SampleID=="MOC",0 ,NA)
data$or <-ifelse(data$SampleID=="1hour",0 ,data$or)
data$or <-ifelse(data$SampleID=="4hour",2 ,data$or)
data$or <-ifelse(data$SampleID=="1day",12 ,data$or)
data$or <-ifelse(data$SampleID=="3day",36 ,data$or)
data$or <-ifelse(data$SampleID=="7day",84 ,data$or)

#MOC only data
moc  <- data[data$SampleID=="MOC",]
#All other time data
time <- data[data$SampleID!="MOC",]

#MOC PLOT
mocplot <-    
    ggplot(moc, aes(x=or, y=mean_size, fill=SampleID)) +
    geom_bar(stat="identity", width =0.2)+
    geom_errorbar(aes(ymin=mean_size, ymax=sd,color=SampleID), width=.2, position=position_dodge(.9)) +
    xlab("") +
    ylab("Mean Distance from PBS")+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_discrete(labels=c("MOC"))+
    scale_fill_manual(values=c("#A88EFF")) +
    scale_color_manual(values=c("#A88EFF")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size=20,face="bold"),panel.spacing.x=unit(2.5, "lines"),axis.title=element_text(size=20, face="bold"))

#All other Times Plot
timeplot <-    
    ggplot(time, aes(x=or, y=mean_size, fill=SampleID)) +
    geom_area(stat = "identity", color= "black",fill="lightgrey", linetype="dashed",alpha=0.2)+
    geom_bar(stat="identity", width =1)+
    geom_errorbar(aes(ymin=mean_size, ymax=sd,color=SampleID), width=.2, position=position_dodge(.9)) +
    xlab("") +
    ylab("")+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_discrete(labels=c("1hour"="1 Hour", "4hour"= "4 Hours",
                                "1day"="1 Day", "3day"="3 Days","7day"="7 Days"))+
    scale_fill_manual(values=c("#6ABD23","#C49A02","#F8766D","#18C59D","#06B9EB")) +
    scale_color_manual(values=c("#6ABD23","#C49A02","#F8766D","#18C59D","#06B9EB")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line.y = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.line.x = element_line(colour = "black"),axis.ticks.x=element_line(colour = "black"),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=unit(2.5, "lines"),
    axis.title=element_text(size=20, face="bold"))

#-------------------
#-------Figure 6F
#-------------------
pdf("Metatrans_InterGroup_Distance_From_PBS.pdf", height = 5, width = 10)
plot_grid(mocplot, timeplot, align = "h", ncol = 2, rel_widths = c(1.7/10, 8.3/10))
dev.off()



#Extract just TAXA Data from Count Table
d2 <- d1[!grepl("UNGROUPED|UNMAPPED",rownames(d1)),]
#keep only data which includes Taxa Data
taxa <- d2[grep("g_|unclassified",rownames(d2)),]
#Create Column with Taxa Names
taxa$gs<- substring(rownames(taxa), 8)
#Create Column with KO NAMES
taxa$ko<- substr(rownames(taxa), 1,6)
#Get DATA just for MOC Innoculum
taxa.moc <- data.frame(MOC=taxa$MOC)
#Keep the rownames
rownames(taxa.moc) <- rownames(taxa)
#Extract the Top 3 TAXA that made up MOC
taxa2 <- taxa[grep("g__Veillonella.s__Veillonella_parvula|g__Streptococcus.s__Streptococcus_mitis_oralis_pneumoniae|g__Prevotella.s__Prevotella_melaninogenica",taxa$gs),]
#Keep only the count data
taxa2 <- taxa2[,1:21]

#Combine count data with MetaDAta
dds.taxa <- DESeqDataSetFromMatrix(countData = taxa2,
                              colData = coldata,
                              design= ~ Time_exp)

#=========================================================
//////////////BAR PLOTS OF SIGNFICANT GENES////////////////
#=========================================================
#Take only count data with Genus Column
taxas <- taxa[,1:22]

#Calculate sum based on Taxa and divide by colsum (relative Abundance)
taxas <-
    taxas %>% 
        group_by(gs) %>% 
        summarise_all(funs(sum)) %>%
        mutate_if(is.numeric, funs(./sum(.))) 

#Melt the Table
as <- data.table(melt(taxas))
#Only keep top 3 TAXA that made up MOC
as <- as[grep("g__Veillonella.s__Veillonella_parvula|g__Streptococcus.s__Streptococcus_mitis_oralis_pneumoniae|g__Prevotella.s__Prevotella_melaninogenica",as$gs),]
#Create a Variable name to match with MetaData
as$SampleID <- as$variable

#Create a Data Frame of MetaData
assay <- as.data.frame(colData(dds.taxa))
assay$SampleID <- assay$Mouse.number

#Merge Assay Data and Meta Data
as2 <- merge(as,assay, all=TRUE)

#Change variable name to Genus
as2$Genus<- as2$gs

#Calculate Median and IQR for each TAXA by the TIME Course
data <- setDT(as2)[,list(value=as.numeric(median(value, na.rm=TRUE)), iqr=as.numeric(quantile(value, probs=.75, na.rm=TRUE))), by=c("Time_exp", "Genus")]

#Remove NEG/n.a.
data <- data[data$Time_exp!="n.a.",]

#Look at the MOC
datamoc <- data[data$Time_exp=="MOC"]
#Create Varibles for Pie Chart
datamoc$ymax = cumsum(datamoc$value)
datamoc$ymin = c(0, head(datamoc$ymax, n=-1))
datamoc <- datamoc[datamoc$value!=0,]
#------------------
#-----------Figure 6A
#------------------
pdf("MetaTrans_MOC_Donut.pdf", height = 5, width = 10)
mocdo <-   
    ggplot(datamoc, aes(fill=Genus, ymax=ymax, ymin=ymin, xmax=8, xmin=3)) +
    geom_rect() +
    coord_polar(theta="y") +
    xlim(c(0, 8)) +
    scale_fill_manual(values=c("#D6504C","#47B665","#FBBE51","#D6D6D6"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.ticks= element_blank(),axis.text=element_blank(),legend.position="none")
dev.off()

#Remove MOC
data2 <- data[data$Time_exp!="MOC"]
data2 <- data2[data2$Time_exp!="PBS"]
#Calculate the Other Taxa rel Abundance
data3 <-
    data2 %>%
    bind_rows(data2 %>% 
            group_by(Time_exp) %>%
            summarise_if(is.numeric, funs(sum))) %>% 
    mutate(Genus = ifelse(is.na(Genus), "Other", Genus))

#Covert Other to 1- and covert all to Percentage
data3$value <- ifelse(data3$Genus=="Other",(1-data3$value)*100,data3$value*100)
data3$or <-ifelse(data3$Time_exp=="1hour",0 ,NA)
data3$or <-ifelse(data3$Time_exp=="4hour",4 ,data3$or)
data3$or <-ifelse(data3$Time_exp=="1day",12 ,data3$or)
data3$or <-ifelse(data3$Time_exp=="3day",36 ,data3$or)
data3$or <-ifelse(data3$Time_exp=="7day",84 ,data3$or)
#------------------
#---------Figure 6G
#------------------
pdf("Metatrans_Relative_Abundance_Median_Stacked_BAL_OTU_over_Time.pdf", height = 7, width = 10)
    ggplot(data3, aes(x=or, y=value)) +
    geom_bar(stat="identity",width=2.5,aes(fill=factor(Genus,
    levels=c("Other","g__Prevotella.s__Prevotella_melaninogenica","g__Streptococcus.s__Streptococcus_mitis_oralis_pneumoniae",
    "g__Veillonella.s__Veillonella_parvula"))))+
    xlab("") +
    ylab("% Relative Abundance")+
    scale_y_continuous(limits=c(0, 105), expand = c(0, 0))+
    scale_x_discrete(labels=c("MOC"="MOC", "1hour"="1Hr", "4hour"= "4Hr",
                                "1day"="1D", "3day"="3D","7day"="7D","PBS"="PBS"))+
    scale_fill_manual(values=c( "#D6D6D6","#D6504C","#47B665","#FBBE51")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank(),
    axis.text.x = element_text(size=20,face="bold"),legend.position = "none",strip.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size=20,face="bold"),panel.spacing.x=unit(2.5, "lines"),axis.title=element_text(size=20, face="bold"))
dev.off()

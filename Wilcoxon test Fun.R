setwd("C:/Users/Kemen Guest User/Documents/Lotus project/Data/FITS2_data")
getwd()
install.packages("dplyr") 

library(tidyr)
library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(funrar)
library(dplyr)
library(lme4)


#####---Table-generation---#####

sampletable <- read.table("Lotara_table_Endo_only.txt" , header = T , check.names = F , stringsAsFactors = F) # mapping file containes informatio about data, lab part
lotus <- subset(sampletable, sampletable$Plant=="Lotus")
lotus_shoot <- subset(lotus, lotus$PlantOrgan=="Shoots")
arabidodpsis <- read.table("arabidopsi_mapfile.txt" , header = T , check.names = F , stringsAsFactors = F)
arabidopsis1 <- subset(arabidodpsis, arabidodpsis$TimePoint=="Fall2017" |  arabidodpsis$TimePoint=="Spring2018" | arabidodpsis$TimePoint=="Fall2018" |  arabidodpsis$TimePoint=="Spring2019")
arabidopsis1$Plant=paste("Arabidopsis")
arabidopsis2 <- subset(arabidopsis1, arabidopsis1$Compartment=="Endo")
arabidopsis2$PlantOrgan=paste("nA")
ara_wc <- arabidopsis2[c(3,1,13,10,4,8,14)]
colnames(ara_wc)[colnames(ara_wc) == "Lib"] <- "Group"
colnames(ara_wc)[colnames(ara_wc) == "TimePoint"] <- "Year"
ara_c = ara_wc
ara_c$Year[ara_wc$Year == "Fall2017" | ara_c$Year == "Spring2018"] <- "2018"
ara_c$Year[ara_wc$Year == "Fall2018" | ara_c$Year == "Spring2019"] <- "2019"

lotara_table <- rbind(ara_c,lotus_shoot)
write.table(lotara_table, file = "Lotara_table_Endo_only.txt",quote = FALSE,sep = "\t",row.names = FALSE)



#####---OTU-Table-generation---#####

Otu = read.table("FITSshared.shared", header = T , check.names = F , stringsAsFactors = F)
taxa = read.table("FITStax.taxonomy", header = T , check.names = F , stringsAsFactors = F)
taxa$Taxonomy = gsub("\\s*\\([^\\)]+\\)","",as.character(taxa$Taxonomy))
otu_taxa <- taxa %>% separate("Taxonomy" , c("Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" ,"Species" , "Strain" ) , ";" , remove = FALSE)
Unclassifiedotus = subset(otu_taxa,otu_taxa$Phylum=="k__Fungi_unclassified" |  otu_taxa$Phylum=="p__unclassified_Fungi") #remove all the otus tat are not classified in phylum level from taxonomy table
Otu1 = Otu[,!(colnames(Otu) %in% Unclassifiedotus$OTU)] #remove remove all the otus tat are not classified in phylum level from otu table
sampletable <- read.table("Lotara_table_Endo_only.txt" , header = T , check.names = F , stringsAsFactors = F) # mapping file containes informatio about data, lab part

Otu2 <- Otu1%>%separate("Group" , c("Loci" , "LibM") , remove = FALSE)
Otu3 = merge(Otu2,sampletable[c("Group","Samplenumber")] , by.x = "LibM" , by.y = "Group")
#Otu2$LibM[!(Otu2$LibM %in% sampleinfo$Lib)]
#Otu3 = merge(sampletable[c("Group","Samplenumber")],Otu2 , by.x = "Group" , by.y = "LibM")
sam_otu <- Otu3[,-c(1:5)]
sam_otu1 = aggregate(. ~  Samplenumber, data = sam_otu, sum)  #some samples are repeted in the lab because of low concentration of DNA so we merged them
#sample
sample = sampletable[!duplicated(sampletable$Samplenumber) , ]  # now remoe duplate sample names
rownames(sample) <- sample$Samplenumber
#adonis befor rarefaction
rownames(sam_otu1) <- sam_otu1$Samplenumber
sam_otu2 = sam_otu1[-1]
#select and save 0 read samples
sam_lessReads = sam_otu2[rowSums(sam_otu2)==0,]
sam_lessReads1 = merge(x =sam_lessReads[,c(1,ncol(sam_lessReads))] , y=sample , by.x="row.names",by.y="row.names")
sam_lessReads1  = sam_lessReads1[,-c(2,2)]
#remove 0 read samples
sam_otu2 = sam_otu2[!rowSums(sam_otu2)==0,] # use this table for Ra log10 and raw reads saving info
#save abundance tabke with samples using sam_otu2
sample1 = sample[(sample$Samplenumber %in% rownames(sam_otu2)),]
rownames(sample1) <- sample1$Samplenumber
####
sam_otu3 = cbind(TotalReads=rowSums(sam_otu2),numOTUs=ncol(sam_otu2),sam_otu2)
cleansam_otu = merge(x = sample1,y = sam_otu3 , by= "row.names")
cleansam_otu = cleansam_otu[-1] ##after removing unclassified otu st phylum level and samples with less than 100 reads
cleansam_otu = cleansam_otu[order(cleansam_otu$TotalReads,decreasing = TRUE),]
taxonomy = otu_taxa[(otu_taxa$OTU %in% colnames(sam_otu3)),]
taxonomy = taxonomy[-ncol(taxonomy)]


#####-WRITE-TABLES-#####

write.table(cleansam_otu,file = "FITS2Otu_ssendo.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(taxonomy,file = "FITS2Taxa_ssendo.txt",quote = FALSE , sep = "\t" , row.names = FALSE)
write.table(sam_lessReads1,file = "FITS2_0readsamples_ssendo.txt",quote = FALSE , sep = "\t" , row.names = FALSE)
write.table(Unclassifiedotus,file = "FITS2Unclassified_ssendo.txt",quote = FALSE , sep = "\t" , row.names = FALSE)

sam_otuRA = make_relative(as.matrix(sam_otu2))          #abundance matrix, with sites in rows and species in columns
t = as.matrix(sam_otu2)
sam_otuRA = data.frame(sam_otuRA)
sam_otuRALog10 = log10(sam_otuRA+1)
sam_otuRALog10 = format(sam_otuRALog10, scientific=F)
sam_otu3RALog10 = cbind(TotalReads=rowSums(sam_otu2),numOTUs=ncol(sam_otu2),sam_otuRALog10)
cleansam_otuRALog10 = merge(x = sample1,y = sam_otu3RALog10 , by= "row.names")
cleansam_otuRALog10 = data.frame(cleansam_otuRALog10[-1]) ##after removing unclassified otu st phylum level and samples with less than 100 reads
cleansam_otuRALog10 = cleansam_otuRALog10[order(cleansam_otuRALog10$TotalReads,decreasing = TRUE, method="auto"),]
write.table(cleansam_otuRALog10,file = "FITS2OtuRALog10_ssendo.txt",quote = FALSE,sep = "\t",row.names = FALSE)



#####---Diversity---#####

cleansam_otu = read.table("FITS2Otu_ssendo.txt", header = T , check.names = F , stringsAsFactors = F)

cleansam_otu1 <- cleansam_otu[order(cleansam_otu$Plant, decreasing=F),]
rarecurve_cleansam <- cleansam_otu1[,-c(1:9)]
cleansam_otu_si <- cleansam_otu1[c(1:9)]
rarecurve(rarecurve_cleansam, step=100, Sample=2, col = c(rep("chartreuse3",127), rep("blue",84)), cex = 0.6,  xlab = "Sample Size", ylab = "Number of OTUs in all samples", label=F)
count(cleansam_otu, cleansam_otu$Plant=="Lotus")

#####RARECURVE-GGPLOT-######


out <- rarecurve(rarecurve_cleansam, step=100, Sample=2, col = c(rep("chartreuse3",127), rep("blue",84)), cex = 0.6,  xlab = "Sample Size", ylab = "Number of OTUs in all samples", label=F)
names(out) <- paste(sample$Plant, sample$Samplenumber, sep = "")


#Let’s plot the results:
par(mfrow=c(1,1))

boxplot(rarefied.si$richness~rarefied.si$Plant, xlab="Plant", ylab="Richness", range=2, col=c("Green","Yellow"), notch=T)
boxplot(rarefied.si$shannon~rarefied.si$Plant, xlab="Plant", ylab="Shannon-Index", range=2, col=c("Green","Yellow"), notch=T)
boxplot(rarefied.si$simpson~rarefied.si$Plant, xlab="Plant", ylab="Simpson-Index", range=2, col=c("Green","Yellow"), notch=T)

boxplot(rarefied.si$richness~rarefied.si$Plant + rarefied.si$Year, xlab="Plant and Year", ylab="Richness", range=2, col=c("Green","Yellow"), notch=T)
boxplot(rarefied.si$shannon~rarefied.si$Plant + rarefied.si$Year , xlab="Plant and Year", ylab="Shannon-Index", range=2, col=c("Green","Yellow"), notch=T)
boxplot(rarefied.si$simpson~rarefied.si$Plant + rarefied.si$Year, xlab="Plant and Year", ylab="Simpson-Index", range=2, col=c("Green","Yellow"), notch=T)

boxplot(rarefied.si$richness~rarefied.si$Plant + rarefied.si$Site, xlab="Plant and Site", ylab="Richness", range=2, col=c("Green","Yellow"))
boxplot(rarefied.si$shannon~rarefied.si$Plant + rarefied.si$Site , xlab="Plant and Site", ylab="Shannon-Index", range=2, col=c("Green","Yellow"))
boxplot(rarefied.si$simpson~rarefied.si$Plant + rarefied.si$Site, xlab="Plant and Site", ylab="Simpson-Index", range=2, col=c("Green","Yellow"))

#####PCOA#####

library(ade4)
?rrarefy
rarefied_rarecurve_cleansam <- rrarefy(rarecurve_cleansam, sample=1000)
vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE)
pca_bray_rarefied.tab<-dudi.pco(vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE), scannf=F)

s.class(pca_bray_rarefied.tab $li, s$Sample, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Bray-Curtis")

#####BRAY-CURTIS#####

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C","#BEAED4" , "#FDC086" , "#FFFF99","#386CB0" , "#F0027F" , "#BF5B17" ,"#999999")

rarefied_rarecurve_cleansam <- rrarefy(rarecurve_cleansam, sample=1000)

richness<-specnumber(rarefied_rarecurve_cleansam, MARGIN = 1)
shannon<-diversity(rarefied_rarecurve_cleansam, index = "shannon", MARGIN = 1, base = exp(1))
simpson<-diversity(rarefied_rarecurve_cleansam, index = "simpson", MARGIN = 1, base = exp(1))
rarefied.si<-cbind(cleansam_otu_si,richness,shannon,simpson)

rarefied.si$Plant = as.factor(rarefied.si$Plant)
rarefied.si$Year = as.factor(rarefied.si$Year)
rarefied.si$Site = as.factor(rarefied.si$Site)
rarefied.si$plantyear =
paste(rarefied.si$Plant, rarefied.si$Year, sep=" ")
rarefied.si$plantyear = as.factor(rarefied.si$plantyear)

vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE)
pca_bray_rarefied<-dudi.pco(vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE), scannf=F)

s.class(pca_bray_rarefied $li, rarefied.si$Plant, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Bray-Curtis")
s.class(pca_bray_rarefied $li, rarefied.si$plantyear, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Bray-Curtis")
s.class(pca_bray_rarefied $li, rarefied.si$Site, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = col2, clabel = 0.5, sub="Bray-Curtis")

#####EUCLIDIAN#####

vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE)
pca_euc_rarefied<-dudi.pco(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), scannf=F)

rarefied.si$Plant = as.factor(rarefied.si$Plant)
rarefied.si$Year = as.factor(rarefied.si$Year)
rarefied.si$Site = as.factor(rarefied.si$Site)

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C","#BEAED4" , "#FDC086" , "#FFFF99","#386CB0" , "#F0027F" , "#BF5B17" ,"#999999")

rarefied_si$Plant
rarefied.si$plantyear = paste(rarefied.si$Plant, rarefied.si$Year, sep="")
rarefied.si$plantyear = as.factor(rarefied.si$plantyear)


s.class(pca_euc_rarefied $li, rarefied.si$Plant, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Euclidean")
s.class(pca_euc_rarefied $li, rarefied.si$plantyear , xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","yellow"), clabel = 0.5, sub="Euclidean")
s.class(pca_euc_rarefied $li, rarefied.si$Site, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = col2, clabel = 0.5, sub="Euclidean")

#####SIMPER#####

inertia(pca_euc_rarefied)

plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "median", members = NULL), labels=rarefied.si$Plant)
plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "complete", members = NULL), labels=rarefied.si$Site)
plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "complete", members = NULL), labels=rarefied.si$Year)



rarefied.si ->rarefied_si
rarefied_si$Plant=as.character(rarefied_si$Plant)
rarefied_si$Year=as.numeric(as.character(rarefied_si$Year))


simper(rarecurve_cleansam, rarefied.si$Plant, permutations = 0, trace = FALSE)
simper(rarecurve_cleansam, rarefied.si$plantyear, permutations = 0, trace = FALSE)
simper(rarecurve_cleansam, rarefied.si$Site, permutations = 0, trace = FALSE)

###WILCOXON-TEST#####
wilcox.test(shannon~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)
wilcox.test(simpson~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)
wilcox.test(richness~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)

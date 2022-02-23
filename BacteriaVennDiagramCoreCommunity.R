####Venn-Diagram####

setwd("C:/Users/Kemen Guest User/Documents/Lotus project/Data/BacV5_data")
getwd()

library(cluster)
library(ape)
require(tidyverse)
install.packages("BiocManager")
BiocManager::install("limma")
library(limma)
install.packages("VennDiagram") 
library("VennDiagram")

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)

#####Read-tables#####


Otu = read.table("BacV5Otu_ssendo.txt", header = T , check.names = F , stringsAsFactors = F)
taxa <- read.table("BacV5Taxa_ssendo.txt" , header = T , check.names = F , stringsAsFactors = F) # mapping file containes informatio about data, lab part
sam_otu1 = Otu
rownames(sam_otu1) = sam_otu1$Samplenumber
sample = sam_otu1[c(1:9)]
sam_otu2 = sam_otu1[-c(1:9)]


#####Split-plants#####

sam_otu_lot = subset(sam_otu1, sam_otu1$Plant=="Lotus")
sam_otu_ara = subset(sam_otu1, sam_otu1$Plant=="Arabidopsis")
sample_lot = sam_otu_lot[c(1:9)]
sample_ara = sam_otu_ara[c(1:9)]
sam_otu_lot2 = sam_otu_lot[-c(1:9)]
sam_otu_ara2 = sam_otu_ara[-c(1:9)]
sam_otu_lotRA <- decostand(sam_otu_lot2, method="total", MARGIN=1)   ###Relative abundance
sam_otu_araRA <- decostand(sam_otu_ara2, method="total", MARGIN=1)   ###Relative abundance
sam_otu_lotRA$Samples = rownames(sam_otu_lotRA)
sam_otu_araRA$Samples = rownames(sam_otu_araRA)


#####--RELATIVE-OCCURRENCE--#####


#####Presence-table#####


threshold = 1 #number more than 0 


#####Transverse#####

sam_otu_ara_tr = t(sam_otu_ara2)
sam_otu_lot_tr = t(sam_otu_lot2)

Lot_Occurrence = rowSums(sam_otu_lot_tr != 0)/ncol(sam_otu_lot_tr) #count non zero columns for each row and devide to number of column
Ara_Occurrence = rowSums(sam_otu_ara_tr != 0)/ncol(sam_otu_ara_tr)


occu_file = cbind(as.data.frame(Lot_Occurrence) , as.data.frame(Ara_Occurrence)) 
occu_file$OTU= paste(rownames(occu_file))
data = merge(occu_file , taxa ,by.x ="OTU" , by.y ="OTU")  #add the information of changing orders name to main data frame


#####-TOP-25-#####


#####Lotus#####


data$Rank_Lotus_RO = rank(-data$Lot_Occurrence, ties.method= "first") 
data_lot25 = subset(data, data$Rank_Lotus_RO<26)
data_lot25$OTUspe = paste(data_lot25$OTU, data_lot25$Species)


#####Arabidopsis#####


data$Rank_Arabidopsis_RO = rank(-data$Ara_Occurrence, ties.method= "first") 
data_ara25 = subset(data, data$Rank_Arabidopsis_RO<26)
data_ara25$OTUspe = paste(data_ara25$OTU, data_ara25$Species)


#####---Matching---#####

RO_match = data_lot25[data_lot25$OTU %in% data_ara25$OTU, ]
RO_nomatch = data_lot25[!data_lot25$OTU %in% data_ara25$OTU, ]

RO_match1 = data_ara25[data_ara25$OTU %in% data_lot25$OTU, ]
RO_nomatch1 = data_ara25[!data_ara25$OTU %in% data_lot25$OTU, ]

RO_match$Rank_Arabidopsis=paste(RO_match1$Rank_Arabidopsis)



#####RELATIVE-ABUNDANCE#####


######--Splitted--######

sam_otu_lot = subset(sam_otu1, sam_otu1$Plant=="Lotus")
sam_otu_ara = subset(sam_otu1, sam_otu1$Plant=="Arabidopsis")
sample_lot = sam_otu_lot[c(1:9)]
sample_ara = sam_otu_ara[c(1:9)]
sam_otu_lot2 = sam_otu_lot[-c(1:9)]
sam_otu_ara2 = sam_otu_ara[-c(1:9)]
sam_otu_lotRA <- decostand(sam_otu_lot2, method="total", MARGIN=1)   ###Relative abundance
sam_otu_araRA <- decostand(sam_otu_ara2, method="total", MARGIN=1)   ###Relative abundance
sam_otu_lotRA$Samples = rownames(sam_otu_lotRA)
sam_otu_araRA$Samples = rownames(sam_otu_araRA)
#melt
sam_otu_araRAMelt = melt(sam_otu_araRA)
sam_otu_lotRAMelt = melt(sam_otu_lotRA)
colnames(sam_otu_lotRAMelt) = c("sample","OTU","RelativeAbundance")
colnames(sam_otu_araRAMelt) = c("sample","OTU","RelativeAbundance")
sam_otu_araRAMelt1 = merge(x =sam_otu_araRAMelt , y = taxa , by.x  = "OTU" , by.y = "OTU")
sam_otu_lotRAMelt1 = merge(x =sam_otu_lotRAMelt , y = taxa , by.x  = "OTU" , by.y = "OTU")


#Lotus
agglo = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$Species, sam_otu_lotRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
agglo = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$OTU, sam_otu_lotRAMelt1$Species, sam_otu_lotRAMelt1$Genus), sum)
colnames(agglo) = c("Species1", "Genus1", "sumAbundance")
colnames(agglo) = c("OTU", "Species1", "Genus1", "sumAbundance")
agglo$Species2 = paste(agglo$Genus1,agglo$Species1, sep="")
agglo$Species2 = paste(agglo$OTU,agglo$Genus1,agglo$Species1, sep=" ")
agglo$Species2[agglo$sumAbundance < 1] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
agglot = subset(agglo, agglo$Species2!="Other")

Lot_comb = merge(data_lot25, agglot,by.x ="OTU" , by.y ="OTU", all=T)
#Lot_comb %>% replace_na(list(Lot_comb$Species = Lot_comb$Species1))
Lot_comb$Species = ifelse(is.na(Lot_comb$Species1) & !is.na(Lot_comb$Species), Lot_comb$Species, Lot_comb$Species1)
Lot_comb$Genus = ifelse(is.na(Lot_comb$Genus1) & !is.na(Lot_comb$Genus), Lot_comb$Genus, Lot_comb$Genus1)
Lot_comb1 = Lot_comb
Lot_comb1[is.na(Lot_comb1)]<-0
Lot_comb1$OTUgenspe= paste(Lot_comb1$OTU, Lot_comb1$Genus, Lot_comb1$Species, sep ="_")


Lot_comb1$Lotus_Core = paste(1)
Lot_comb1$Lotus_Core[Lot_comb1$Lot_Occurrence < 0] <- 1
Lot_comb1$Lotus_Core[Lot_comb1$Lot_Occurrence==0] <- 0

Lot_comb1$Lotus_Abundance = paste(1)
Lot_comb1$Lotus_Abundance[Lot_comb1$sumAbundance==0] <- 0

Lot_comb1= Lot_comb1[c(19,20,21)]


#Arabidopsis
aggar = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
aggar = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$OTU, sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum)
colnames(aggar) = c("Species1", "Genus1", "sumAbundance")
colnames(aggar) = c("OTU","Species1", "Genus1", "sumAbundance")
aggar$Species2 = paste(aggar$Genus1,aggar$Species1, sep="")
aggar$Species2 = paste(aggar$OTU,aggar$Genus1,aggar$Species1, sep=" ")
#agg$Species2 <- str_replace(agg$Species2, "_unclassified", "Other")
aggar$Species2[aggar$sumAbundance < 1.7] <- "Other"##change the order of otus that have less than a threshold to "Other"
aggara = subset(aggar, aggar$Species2!="Other")
Ara_comb = merge(data_ara25 , aggara ,by.x ="OTU" , by.y ="OTU", all=T)

Ara_comb$Species = ifelse(is.na(Ara_comb$Species1) & !is.na(Ara_comb$Species), Ara_comb$Species, Ara_comb$Species1)
Ara_comb$Genus = ifelse(is.na(Ara_comb$Genus1) & !is.na(Ara_comb$Genus), Ara_comb$Genus, Ara_comb$Genus1)
Ara_comb1 = Ara_comb
Ara_comb1[is.na(Ara_comb1)]<-0
Ara_comb1$OTUgenspe= paste(Ara_comb1$OTU, Ara_comb1$Genus, Ara_comb1$Species, sep ="_")


Ara_comb1$Arabidopsis_Core = paste(1)
Ara_comb1$Arabidopsis_Core[Ara_comb1$Ara_Occurrence < 0] <- 1
Ara_comb1$Arabidopsis_Core[Ara_comb1$Ara_Occurrence==0] <- 0

Ara_comb1$Arabidopsis_Abundance = paste(1)
Ara_comb1$Arabidopsis_$Abundance[Ara_comb1$sumAbundance==0] <- 0

Ara_comb1= Ara_comb1[c(20,21,22)]

completetable = merge(Ara_comb1 , Lot_comb1 ,by.x ="OTUgenspe" , by.y ="OTUgenspe", all=T)
completetable[is.na(completetable)]<-0

completetable1 = completetable[c(2,3,4,5)]
completetable1$Arabidopsis_Core = as.numeric(completetable1$Arabidopsis_Core)
completetable1$Arabidopsis_Abundance = as.numeric(completetable1$Arabidopsis_Abundance)
completetable1$Lotus_Core = as.numeric(completetable1$Lotus_Core)
completetable1$Lotus_Abundance = as.numeric(completetable1$Lotus_Abundance)
completetable2 = vennCounts(completetable1)
vennDiagram(completetable2, names = c("Arabidopsis Abundance", "Arabidopsis Core", "Lotus Abundance", "Lotus Core"), cex = 1, counts.col = "red", circle.col = col)
col = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#416B4C","#CE5A17","#776804","#2B2B83","#416B4C","#CE5A17","#776804","#2B2B83","#EC2B80" ,"#EDB8B6" ,"#AB56E7","#999999")                           
?par

#Tubingen data analyses 
packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)
getwd()
install.packages("ape")
library(cluster)
library(ape)
library(tidyverse)
library(tidyr)

#path this project : ""C:/Users/Kemen Guest User/Documents/Lotus project/Data/BacV5_data""
#load packages
setwd("C:/Users/Kemen Guest User/Documents/Lotus project/Data/BacV5_data")
getwd()
#Bacteria
Otu = read.table("BacV5Otu_ssendo.txt", header = T , check.names = F , stringsAsFactors = F)
taxa <- read.table("BacV5Taxa_ssendo.txt" , header = T , check.names = F , stringsAsFactors = F) # mapping file containes informatio about data, lab part
sam_otu1 = Otu
rownames(sam_otu1) = sam_otu1$Samplenumber
sample = sam_otu1[c(1:9)]
sam_otu2 = sam_otu1[-c(1:9)]
sam_otuRA <- decostand(sam_otu2, method="total", MARGIN=1)   ###Relative abundance
sam_otuRA$Samples = rownames(sam_otuRA)
#melt
sam_otuRAMelt = melt(sam_otuRA)
colnames(sam_otuRAMelt) = c("sample","OTU","RelativeAbundance")
sam_otuRAMelt1 = merge(x =sam_otuRAMelt , y = taxa , by.x  = "OTU" , by.y = "OTU")


#####---Phylum---#####


agg = aggregate(sam_otuRAMelt1[,3], list(sam_otuRAMelt1$Phylum), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg) = c("Phylum1","sumAbundance")
agg$Phylum2 = agg$Phylum1  # copy of this col because i want to have main col and change the order of this column to "other" category

agg$Phylum2 = as.character(agg$Phylum2)
agg$Phylum2[agg$sumAbundance < 0.5] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
data = merge(sam_otuRAMelt1 , agg ,by.x ="Phylum" , by.y ="Phylum1")  #add the information of changing orders name to main data frame
data$Phylum2 <- reorder(data$Phylum2 , data$sumAbundance)
data1 = merge(x =data , y = sample , by.x  = "sample" , by.y = "Samplenumber")

col = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#416B4C","#CE5A17","#776804","#2B2B83","#416B4C","#CE5A17","#776804","#2B2B83","#EC2B80" ,"#EDB8B6" ,"#AB56E7","#999999")

data1$Year = as.factor(data1$Year)
data1$Plant = as.factor(data1$Plant)
data1$Site = as.factor(data1$Site)

Ra = "RelativeAbundance"
Phylum = "Phylum2"
variable = c("Plant","Year","Site")
for (i in 1:3){
  var = variable[i]
  a <- ggplot()  + geom_bar(data=data1, aes_string(x=var, y=Ra, fill = Phylum), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  +
    scale_fill_manual(values=rev(col)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_bw() + theme(text = element_text(size=30 , colour = "black"),
                       axis.text = element_text(size = 30 , colour = "black" ),
                       axis.title = element_text(size = 30 , colour = "black" ),
                       strip.text = element_text(size = 30 , colour = "black" ),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black") ,
                       legend.title = element_text(color = "black", size = 32),
                       legend.text = element_text(color = "black", size = 30) ,
                       strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Phylum") +
    ggsave(paste("BacV5Phylum",variable[i],".png" , sep = "") , a , width = 25 , height = 8 , units = "in")
  ggsave(paste("BacV5Phylum",variable[1],".png" , sep = ""), a , width = 45 , height = 8 , units = "in")
  
}


#####---Order---#####


agg = aggregate(sam_otuRAMelt1[,3], list(sam_otuRAMelt1$Order), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg) = c("Order1","sumAbundance")
agg$Order2 = agg$Order1  # copy of this col because i want to have main col and change the order of this column to "other" category
agg$Order2 = as.character(agg$Order2)
sort(agg$sumAbundance)
agg$Order2[agg$sumAbundance < 	5.8] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
data = merge(sam_otuRAMelt1 , agg ,by.x ="Order" , by.y ="Order1")  #add the information of changing orders name to main data frame
data$Order2 <- reorder(data$Order2 , data$sumAbundance)
data1 = merge(x =data , y = sample , by.x  = "sample" , by.y = "Samplenumber")


data1$Year = as.factor(data1$Year)
data1$Plant = as.factor(data1$Plant)
data1$Site = as.factor(data1$Site)
data1$Order = as.factor(data1$Order)

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#999999")

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C","#BEAED4" , "#FDC086" , "#FFFF99" ,"#999999")


Ra = "RelativeAbundance"
Order2 = "Order2"
variable = c("Plant","Year","Site")


for (i in 1:3){
  var = variable[i]
  p <- ggplot()  + geom_bar(data=data1, aes_string(x=var, y=Ra, fill = Order2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col2)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Order") +
    #ggsave("Results/FungiorderAll1.pdf" , p , width = 15 , height = 10 , units = "in")
    ggsave(paste("BacV5Order",variable[i],".png" , sep = "") , p , width = 15, height = 10, units = "in")
  ggsave(paste("BacV5Order",variable[1],".png" , sep = "") , p , width = 45 , height = 10 , units = "in")
}  


#####---Genus---#####


agg = aggregate(sam_otuRAMelt1[,3], list(sam_otuRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg) = c("Genus1","sumAbundance")
agg$Genus2 = agg$Genus1  # copy of this col because i want to have main col and change the order of this column to "other" category
agg$Genus2 = as.character(agg$Genus2)
sort(agg$sumAbundance)
agg$Genus2[agg$sumAbundance < 5] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
data = merge(sam_otuRAMelt1 , agg ,by.x ="Genus" , by.y ="Genus1")  #add the information of changing orders name to main data frame
data$Genus2 <- reorder(data$Genus2 , data$sumAbundance)
data1 = merge(x =data , y = sample , by.x  = "sample" , by.y = "Samplenumber")


col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")


variable = c("Plant","Year","Site")

data1$Year = as.factor(data1$Year)
data1$Plant = as.factor(data1$Plant)
data1$Site = as.factor(data1$Site)
Genus2 = "Genus2"
Ra = "RelativeAbundance"

for (i in 1:3){
  var = variable[i]
  z <- ggplot()  + geom_bar(data=data1, aes_string(x=var, y=Ra, fill = Genus2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col2)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Genus") +
    ggsave(paste("BacV5Genus",variable[i],".png" , sep = "") , z , width = 15, height = 10, units = "in")
  ggsave(paste("BacV5Genus",variable[1],".png" , sep = "") , z , width = 45 , height = 10 , units = "in")
}  


#####---Species---#####


agg = aggregate(sam_otuRAMelt1[,3], list(sam_otuRAMelt1$Species, sam_otuRAMelt1$Genus, sam_otuRAMelt1$OTU), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg) = c("Species1", "Genus1", "OTU", "sumAbundance")
agg$Species2 = paste(agg$OTU,agg$Genus1,agg$Species1, sep=" ")
agg$Species2 = as.character(agg$Species2)
agg$Species1 = as.character(agg$Species1)
sort(agg$sumAbundance)
agg=agg[c(1,4,5)]
#agg$Species2 <- str_replace(agg$Species2, "_unclassified", "Other")
?str_replace
#agg$Species2 <- gsub('_unclassified', "Other", agg$Species1)
agg$Species2[agg$sumAbundance < 2.2] <- "Other"##change the order of otus that have less than a threshold to "Other"
data = merge(sam_otuRAMelt1 , agg , by.x ="Species" , by.y ="Species1")  #add the information of changing orders name to main data frame
data$Species2 <- reorder(data$Species2 , data$sumAbundance)
data1 = merge(x =data , y = sample , by.x  = "sample" , by.y = "Samplenumber")



col3 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C","#BEAED4" , "#FDC086" , "#FFFF99","#386CB0" , "#F0027F" , "#BF5B17" ,"#999999")

#####--Merged--#####


variable = c("Plant","Year","Site")
data1$Year = as.factor(data1$Year)
data1$Plant = as.factor(data1$Plant)
data1$Site = as.factor(data1$Site)
Species2="Species2"
Ra = "RelativeAbundance"

for (i in 1:3){
  var = variable[i]
  y <- ggplot()  + geom_bar(data=data1, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    ggsave(paste("BacV5Species",variable[i],".png" , sep = "") , y , width = 25, height = 12, units = "in")
  ggsave(paste("BacV5Species",variable[1],".png" , sep = "") , y , width = 45 , height = 12 , units = "in")
}  


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
agg = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$Species, sam_otu_lotRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
agg = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$OTU, sam_otu_lotRAMelt1$Species, sam_otu_lotRAMelt1$Genus), sum)
colnames(agg) = c("Species1", "Genus1", "sumAbundance")
colnames(agg) = c("OTU", "Species1", "Genus1", "sumAbundance")
agg$Species2 = paste(agg$Genus1,agg$Species1, sep="")
agg$Species2 = paste(agg$OTU,agg$Genus1,agg$Species1, sep=" ")
agg$Species2[agg$sumAbundance < 1] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
data_lo = merge(sam_otu_lotRAMelt1 , agg ,by.x ="OTU" , by.y ="OTU")  #add the information of changing orders name to main data frame
data_lo$Species2 <- reorder(data_lo$Species2 , data_lo$sumAbundance)
data_lot = merge(x =data_lo , y = sample_lot , by.x  = "sample" , by.y = "Samplenumber")
data_lot$Year = as.factor(data_lot$Year)
Species2="Species2"
Ra = "RelativeAbundance"
variable = c("Plant","Year","Site")


for (i in 1:3){
  var = variable[i]
  t <- ggplot()  + geom_bar(data=data_lot, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    #ggsave("Results/FungiorderAll1.pdf" , p , width = 15 , height = 10 , units = "in")
    ggsave(paste("LotusBacteriaSpecies",variable[i],".png" , sep = "") , t , width = 25, height = 12, units = "in")
  ggsave(paste("LotusBacteriaSpecies",variable[1],".png" , sep = "") , t , width = 45 , height = 12 , units = "in")
}  

#Arabidopsis
agg = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
agg = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$OTU, sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum)
colnames(agg) = c("Species1", "Genus1", "sumAbundance")
colnames(agg) = c("OTU","Species1", "Genus1", "sumAbundance")
agg$Species2 = paste(agg$Genus1,agg$Species1, sep="")
agg$Species2 = paste(agg$OTU,agg$Genus1,agg$Species1, sep=" ")
#agg$Species2 <- str_replace(agg$Species2, "_unclassified", "Other")
?str_replace
agg$Species2 <- gsub('_unclassified', " sp.", agg$Species2)
agg$Species2 <- gsub('g__', "", agg$Species2)
agg$Species2 <- gsub('f__', "", agg$Species2)
agg$Species2 <- gsub('s__', " ", agg$Species2)
agg$Species2[agg$sumAbundance < 1.7] <- "Other"##change the order of otus that have less than a threshold to "Other"
data_ar = merge(sam_otu_araRAMelt1 , agg ,by.x ="OTU" , by.y ="OTU")  #add the information of changing orders name to main data frame
data_ar$Species2 <- reorder(data_ar$Species2 , data_ar$sumAbundance)
data_ara = merge(x =data_ar , y = sample_ara , by.x  = "sample" , by.y = "Samplenumber")
data_ara$Year = as.factor(data_ara$Year)
Species2="Species2"
Ra = "RelativeAbundance"
variable = c("Plant","Year","Site")

for (i in 1:3){
  var = variable[i]
  v <- ggplot()  + geom_bar(data=data_ara, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    #ggsave("Results/FungiorderAll1.pdf" , p , width = 15 , height = 10 , units = "in")
    ggsave(paste("ArabidopsisBacteriaSpecies",variable[i],".png" , sep = "") , v , width = 25, height = 12, units = "in")
  ggsave(paste("ArabidopsisBacteriaSpecies",variable[1],".png" , sep = "") , v , width = 45 , height = 12 , units = "in")
}  


#####---Sort-by-RA---#####


#####Lotus#####

agg_lo = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$Species, sam_otu_lotRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg_lo) = c("Species1", "Genus1", "sumAbundance")
agg_lo$Species2 = paste(agg_lo$Genus1,agg_lo$Species1, sep="")

agg_lo$Rank_Lotus = rank(-agg_lo$sumAbundance, ties.method= "min") 
agg_lo50 = subset(agg_lo, agg_lo$Rank_Lotus<51)
agg_lo50$spra = paste(agg_lo50$Rank_Lotus, agg_lo50$Species2)


#####Arabidopsis#####

agg_ar = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum) ## sum of relative abundance of different orders in all samples colnames(agg) = c("order1" , "sumAbundance")
colnames(agg_ar) = c("Species1", "Genus1", "sumAbundance")
agg_ar$Species2 = paste(agg_ar$Genus1,agg_ar$Species1, sep="")

agg_ar$Rank_Arabidopsis = rank(-agg_ar$sumAbundance, ties.method= "min") 
agg_ar50 = subset(agg_ar, agg_ar$Rank_Arabidopsis<51)
agg_ar50$spra = paste(agg_ar50$Rank_Arabidopsis, agg_ar50$Species2)


#####---Matching---#####

agg_match = agg_lo50[agg_lo50$Species2 %in% agg_ar50$Species2, ]
agg_nomatch = agg_lo50[!agg_lo50$Species2 %in% agg_ar50$Species2, ]

agg_match1 = agg_ar50[agg_ar50$Species2 %in% agg_lo50$Species2, ]
agg_nomatch1 = agg_ar50[!agg_ar50$Species2 %in% agg_lo50$Species2, ]

agg_match$Rank_Arabidopsis=paste(agg_match1$Rank_Arabidopsis)
agg_nomatch$Rank_Arabidopsis=paste(agg_nomatch1$Rank_Arabidopsis)


#####-Write-tables#####

write.table(agg_match,file = "BacteriaSpeciesOverlap.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(agg_nomatch1,file = "BacteriaNoMatchArabidopsis.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(agg_nomatch,file = "BacteriaNoMatchLotus.txt",quote = FALSE,sep = "\t",row.names = FALSE)


#####---Plotting---#####


#####Matching-results#####

data_m = merge(sam_otuRAMelt1 , agg_match ,by.x ="Species" , by.y ="Species1")  #add the information of changing orders name to main data frame
data_m$Species2 <- reorder(data_m$Species2 , data_m$sumAbundance)
data_ma = merge(x =data_m , y = sample , by.x  = "sample" , by.y = "Samplenumber")

col3 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")

variable = c("Plant","Year","Site")
data_ma$Year = as.factor(data_ma$Year)
data_ma$Plant = as.factor(data_ma$Plant)
data_ma$Site = as.factor(data_ma$Site)
Species2="Species2"
Ra = "RelativeAbundance"

for (i in 1:3){
  var = variable[i]
  c <- ggplot()  + geom_bar(data=data_ma, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    ggsave(paste("MatchingBacteriaSpecies",variable[i],".png" , sep = "") , c , width = 25, height = 12, units = "in")
  ggsave(paste("MatchingBacteriaSpecies",variable[1],".png" , sep = "") , c , width = 45 , height = 12 , units = "in")
}  


#####Non-matching-results#####

#####Lotus#####

data_lotn = merge(sam_otu_lotRAMelt1 , agg_nomatch ,by.x ="Species" , by.y ="Species1")  #add the information of changing orders name to main data frame
data_lotn$Species2 <- reorder(data_lotn$Species2 , data_lotn$sumAbundance)
data_lotnm = merge(x =data_lotn , y = sample , by.x  = "sample" , by.y = "Samplenumber")

col3 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")

variable = c("Plant","Year","Site")
data_lotnm$Year = as.factor(data_lotnm$Year)
Species2="Species2"
Ra = "RelativeAbundance"

for (i in 1:3){
  var = variable[i]
  b <- ggplot()  + geom_bar(data=data_lotnm, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    ggsave(paste("NonMatchingLotusBacteriaSpecies",variable[i],".png" , sep = "") , b , width = 25, height = 12, units = "in")
  ggsave(paste("NonMatchingLotusBacteriaSpecies",variable[1],".png" , sep = "") , b , width = 45 , height = 12 , units = "in")
}  


#####Arabidopsis#####

data_aran = merge(sam_otu_araRAMelt1 , agg_nomatch1 ,by.x ="Species" , by.y ="Species1")  #add the information of changing orders name to main data frame
data_aran$Species2 <- reorder(data_aran$Species2 , data_aran$sumAbundance)
data_aranm = merge(x =data_lotn , y = sample , by.x  = "sample" , by.y = "Samplenumber")

col3 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")

variable = c("Plant","Year","Site")
data_aranm$Year = as.factor(data_aranm$Year)
Species2="Species2"
Ra = "RelativeAbundance"

for (i in 1:3){
  var = variable[i]
  d <- ggplot()  + geom_bar(data=data_lotnm, aes_string(x=var, y=Ra, fill = Species2), stat="identity", position="fill") +
    facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    guides(fill=guide_legend(ncol=1))  + 
    theme_classic() + theme(text = element_text(size=30 , colour = "black"),
                            axis.text = element_text(size = 30 , colour = "black" ),
                            axis.title = element_text(size = 30 , colour = "black" ),
                            strip.text = element_text(size = 30 , colour = "black" ),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), 
                            axis.line = element_line(colour = "black") ,
                            legend.title = element_text(color = "black", size = 30),
                            legend.text = element_text(color = "black", size = 30) ,
                            strip.background = element_rect(colour = "black", fill = "white")) +
    ylab("Relative abundance") + 
    xlab(var) + labs(fill = "Species") +
    ggsave(paste("NonMatchingArabidopsisBacteriaSpecies",variable[i],".png" , sep = "") , d , width = 25, height = 12, units = "in")
  ggsave(paste("NonMatchingArabidopsisBacteriaSpecies",variable[1],".png" , sep = "") , d , width = 45 , height = 12 , units = "in")
}  

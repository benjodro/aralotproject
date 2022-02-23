setwd("C:/Users/Kemen Guest User/Documents/Lotus project/Data/BacV5_data")
getwd()

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)
library(lme4)


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
lotara_RO = as.data.frame(t(occu_file))
core = occu_file[rowSums(occu_file >= threshold) >= experiment , ] #occurence in more than threshold  in each year
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


#####--RELATIVE-ABUNDANCE--######

######Splitted######


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

agglot = aggregate(sam_otu_lotRAMelt1[,3], list(sam_otu_lotRAMelt1$OTU, sam_otu_lotRAMelt1$Genus, sam_otu_lotRAMelt1$Species), sum)
colnames(agglot) = c("OTU","Species1","Genus1","sumAbundance")
agglot$Species2 = paste(agglot$OTU,agglot$Species1,agglot$Species1, sep=" ")
agglot$Rank_Lotus_RA = rank(-agglot$sumAbundance, ties.method= "first") 
agglot25 = subset(agglot, agglot$Rank_Lotus_RA<26)


overlaplot_merge = merge(x=data_lot25 , y = agglot25 , by.x  = "OTU" , by.y = "OTU")
overlaplot_merge1 = overlaplot_merge[c(1,2,3,17,4,18,14,19,13)]

#Arabidopsis

aggara = aggregate(sam_otu_araRAMelt1[,3], list(sam_otu_araRAMelt1$OTU, sam_otu_araRAMelt1$Species, sam_otu_araRAMelt1$Genus), sum)
colnames(aggara) = c("OTU","Species1", "Genus1", "sumAbundance")
aggara$Species2 = paste(aggara$OTU,aggara$Genus1,aggara$Species1, sep=" ")
aggara$Rank_Arabidopsis_RA = rank(-aggara$sumAbundance, ties.method= "first") 
aggara25 = subset(aggara, aggara$Rank_Arabidopsis_RA<26)


overlapara_merge = merge(x = data_ara25 , y = aggara25 , by.x  = "OTU" , by.y = "OTU")
overlapara_merge1 = overlapara_merge[c(1,2,3,18,4,19,15,20,14)]


#####---Write-tables---#####


write.table(overlaplot_merge1,file = "BacteriaLotusSpeciesTOP25_RA_RO.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(overlapara_merge1,file = "BacteriaArabidopsisTOP25_RA_RO.txt",quote = FALSE,sep = "\t",row.names = FALSE)



#####Plotting#####



variable = c("Lot_Occurrence","Ara_Occurrence")
var1 = c("Rank_Lotus","Rank_Arabidopsis")
overlaplot_merge1$Rank_Lotus_RO = as.factor(overlaplot_merge1$Rank_Lotus_RO)
overlapara_merge1$Rank_Arabidopsis_RO = as.factor(overlapara_merge1$Rank_Arabidopsis_RO)
Species2="Species2"
OTUspe="OTUspe"



col3 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C" , "#BEAED4" , "#FDC086" , "#FFFF99" , "#386CB0" , "#F0027F" , "#BF5B17" , "#EDB8B6" , "#1B9E77" , "#D95F02" ,
         "#7570B3" , "#E7298A" , "#66A61E" , "#E6AB02" , "#A6761D" , "#666666" , "#A6CEE3" , "#1F78B4","#B2DF8A" , "#33A02C" ,
         "#FB9A99" , "#E31A1C" , "#FDBF6F" , "#FF7F00" , "#CAB2D6" , "#6A3D9A" , "#AB56E7" , "#B15928" , "#FBB4AE" , "#B3CDE3" , 
         "#CCEBC5" , "#DECBE4" , "#FED9A6" ,  "#FFFFCC" , "#E5D8BD" , "#FDDAEC" ,"#F2F2F2" , "#B3E2CD" , "#FDCDAC" , "#CBD5E8" ,
         "#F4CAE4" , "#E6F5C9" , "#FFF2AE" ,"#CCCCCC" ,"#E41A1C", "#377EB8" ,"#4DAF4A" ,"#999999")

scale = function(x) sprintf("%.3f",x)


#####--Overlap-Lotus--#####

for (i in 1:2){
  var = variable[1]
  b <- ggplot()  + geom_bar(data=overlaplot_merge1, aes_string(x=Species2, y=var, fill = Species2), stat="identity") +
    #facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    scale_y_continuous(labels = scale,limits = c(0,1)) + 
    scale_x_discrete(limit = overlaplot_merge1$Species2, labels = overlaplot_merge1$Rank_Lotus_RO) + 
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
    ylab("Relative Occurrence") + 
    xlab("Rank Lotus") +
    ggsave(paste("Overlapping_RO_RA_BacteriaSpecies",variable[1],".png" , sep = "") , b , width = 25, height = 15, units = "in")
  #ggsave(paste("MatchingROBacteriaSpecies",variable[1],".png" , sep = "") , a , width = 45 , height = 12 , units = "in")
}

#####--Overlap-Arabidopsis--#####

for (i in 1:2){
  var = variable[2]
  b <- ggplot()  + geom_bar(data=overlapara_merge1, aes_string(x=Species2, y=var, fill = Species2), stat="identity") +
    #facet_wrap(.~ Plant, scales = 'free_x' , nrow = 1 , strip.position="bottom")  + 
    scale_fill_manual(values=rev(col3)) +
    scale_y_continuous(labels = scale,limits = c(0,1)) + 
    scale_x_discrete(limit = overlapara_merge1$Species2, labels = overlapara_merge1$Rank_Arabidopsis_RO) + 
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
    ylab("Relative Occurrence") + 
    xlab("Rank Arabidopsis") +
    ggsave(paste("Overlapping_RO_RA_BacteriaSpecies",variable[2],".png" , sep = "") , b , width = 25, height = 15, units = "in")
  #ggsave(paste("MatchingROBacteriaSpecies",variable[1],".png" , sep = "") , a , width = 45 , height = 12 , units = "in")
}

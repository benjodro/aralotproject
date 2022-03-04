setwd("/Users/mahmoudi/Documents/Doc/2020/BengtProject/qpcr/")


library(cluster)
library(ape)
require(tidyverse)

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)

rawcq <- read.csv("Cq values.csv" , header = T , check.names = F , stringsAsFactors = F, sep="\t") # mapping file containes informatio about data, lab part
rawcq = rawcq[-c(17,25,26),]
rawcq1  = separate(rawcq, `PCRproduct-treatment`,into = c("group","Sample") , remove = F , sep = "-" )
rawcq1$Sample = factor(rawcq1$Sample , levels = c("NTC" , "C2" , "S2"))
racont = subset(rawcq1 , rawcq1$Sample == "C2")
racont$Cq_Ratio = racont$`Cq mean`/19.76710

rasam = subset(rawcq1 , rawcq1$Sample == "S2")
rasam$Cq_Ratio = rasam$`Cq mean`/19.25116
ntc = subset(rawcq1 , rawcq1$Sample == "NTC")
ntc$Cq_Ratio = 0


all = rbind(ntc,rasam,racont )
all$PCRproduct_treatment = all$`PCRproduct-treatment` 
#allun = all[!duplicated(all[ , c("PCRproduct_treatment")]), ]

all$ID1 = paste(all$PCRproduct_treatment,1:nrow(all) , sep = "_")
rownames(all) = all$ID1
sdall = all[c(7,4)]


library(dplyr)
ag = sdall %>% group_by(PCRproduct_treatment) %>% summarise_each(funs(mean, sd))

#colnames(ag) = c("PCRproduct_treatment","Cq.mean" , "Cq.sd") 
all1 = merge(all ,ag , by.x = "PCRproduct_treatment"  , by.y = "PCRproduct_treatment" ,)
#all1[is.na(all1)] <- 0
all2 = all1[!duplicated(all1[ , c("PCRproduct_treatment")]), ]

all2$PCRproduct_treatment = factor(all2$PCRproduct_treatment , levels = c("LotEF1a-NTC" , "LotEF1a-C2" , "LotEF1a-S2" , 
                                                                          "Brev-NTC"  ,"Brev-C2"  , "Brev-S2" ,
                                                                          "Dio-NTC" , "Dio-C2" , "Dio-S2"  , 
                                                                          "Cysto-NTC" , "Cysto-C2"  , "Cysto-S2" ))

p <- ggplot(all2, aes(x=PCRproduct_treatment, y=Cq_Ratio , fill = Sample)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Cq_Ratio-sd, ymax=Cq_Ratio+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust=1) ,text = element_text(size=15)) + xlab("Treatment")

ggsave("cq_ratioplot.png" , p  , width = 7 , height =5 )











#rawcq$Cq[rawcq$Cq=="NaN"] <- 0

rawcq$cqmean = mean(rawcq$Cq[rawcq$`PCRproduct-treatment`=="LotEF1a-C2"])
rawcq$cqratiocontrol = rawcq$"Cq mean"/(rawcq$"Cq mean"[1])
rawcq$cqratiotreated = rawcq$"Cq mean"/(rawcq$"Cq mean"[4])

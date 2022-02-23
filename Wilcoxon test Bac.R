setwd("C:/Users/Kemen Guest User/Documents/Lotus project/Data/BacV5_data")
getwd()

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" , "gridExtra" , "funrar" , "ade4" , 
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)




#####---Diversity---#####

cleansam_otu = read.table("BacV5Otu_ssendo.txt", header = T , check.names = F , stringsAsFactors = F)

cleansam_otu1 <- cleansam_otu[order(cleansam_otu$Plant, decreasing=F),]
rarecurve_cleansam <- cleansam_otu1[,-c(1:9)]
cleansam_otu_si <- cleansam_otu1[c(1:9)]
rarecurve(rarecurve_cleansam, step=100, Sample=2, col = c(rep("chartreuse3",124), rep("blue",83)), cex = 0.6,  xlab = "Sample Size", ylab = "Number of OTUs in all samples", label=F)
sort(cleansam_otu1$Plant)

#####CHECKING-DATA-COUNT

cleansam_otu$Plant = as.factor(cleansam_otu$Plant)
cleansam_otu %<% count(Plant)
count(cleansam_otu, cleansam_otu$Plant=="Arabidopsis")

#####RARECURVE-GGPLOT-######

sampletable <- read.table("Lotara_table.txt" , header = T , check.names = F , stringsAsFactors = F) # mapping file containes informatio about data, lab part
sample = sampletable[!duplicated(sampletable$Samplenumber) , ]  # now remoe duplate sample names
zeroread = read.table("BacV5_0readsamples_ss.txt", header = T , check.names = F , stringsAsFactors = F)
sample = sample[,!(colnames(sample) %in% zeroread$Samplenumber)]


out <- rarecurve(rarecurve_cleansam, step=100, Sample=2, col = c(rep("chartreuse3",251), rep("blue",83)), cex = 0.6,  xlab = "Sample Size", ylab = "Number of OTUs in all samples", label=F)
names(out) <- paste(sample$Plant, sample$Samplenumber, sep = "")

?rbind
?mapply
?ggplot

# Coerce data into "long" form.
protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x, stringsAsFactors = default.stringsAsFactors())
  colnames(mydf) <- "value"
  mydf$species <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = FALSE)

xy <- do.call(rbind, protox)
rownames(xy) <- NULL  # pretty

# Plot.
ggplot(xy, aes(x = subsample, y = value, color = species)) +
  theme_bw() +
  ylab("Number of OTU") +
  xlab("Sample size") +
  scale_color_discrete(guide = F) +  # turn legend on or off
  geom_line()

#####-BOXPLOT-GENERATION-#####


rarefied_rarecurve_cleansam <- rrarefy(rarecurve_cleansam, sample=1000)

richness<-specnumber(rarefied_rarecurve_cleansam, MARGIN = 1)
shannon<-diversity(rarefied_rarecurve_cleansam, index = "shannon", MARGIN = 1, base = exp(1))
simpson<-diversity(rarefied_rarecurve_cleansam, index = "simpson", MARGIN = 1, base = exp(1))
rarefied.si<-cbind(cleansam_otu_si,richness,shannon,simpson)

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

col2 = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#7FC97F","#CE5A17","#776804","#2B2B83",
         "#416B4C","#BEAED4" , "#FDC086" , "#FFFF99","#386CB0" , "#F0027F" , "#BF5B17" ,"#999999")

rarefied.si$plantyear = paste(rarefied.si$Plant, rarefied.si$Year, sep=" ")
rarefied.si$plantyear = as.factor(rarefied.si$plantyear)
rarefied.si$Plant = as.factor(rarefied.si$Plant)
rarefied.si$Year = as.factor(rarefied.si$Year)
rarefied.si$Site = as.factor(rarefied.si$Site)

#####BRAY-CURTIS#####

vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE)
pca_bray_rarefied<-dudi.pco(vegdist(rarefied_rarecurve_cleansam, method="bray", binary=FALSE), scannf=F)

s.class(pca_bray_rarefied $li, rarefied.si$Plant, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Bray-Curtis")
s.class(pca_bray_rarefied $li, rarefied.si$plantyear, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Bray-Curtis")
s.class(pca_bray_rarefied $li, rarefied.si$Site, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = col2, clabel = 0.5, sub="Bray-Curtis")

#####EUCLIDIAN#####

vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE)
pca_euc_rarefied<-dudi.pco(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), scannf=F)


s.class(pca_euc_rarefied $li, rarefied.si$Plant, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","brown"), clabel = 0.5, sub="Euclidean")
s.class(pca_euc_rarefied $li, rarefied.si$plantyear , xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = c("chartreuse3","blue","red","yellow"), clabel = 0.5, sub="Euclidean")
s.class(pca_euc_rarefied $li, rarefied.si$Site, xax=1, yax=2, cpoint=3, grid=F, addaxes=T, cellipse=1, cstar=0, axesell=0, col = col2, clabel = 0.5, sub="Euclidean")


#####LMER??#####

lmer (rarefied_si$Year + rarefied_si$Plant, sub="Euclidean")


#####----H-CLUSTER-DENDROGRAM-----#####

inertia(pca_euc_rarefied)
inertia(pca_bray_rarefied)

plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "complete", members = NULL), labels=rarefied.si$Plant)
plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "complete", members = NULL), labels=rarefied.si$Site)
plot(hclust(vegdist(rarefied_rarecurve_cleansam, method="euclidean", binary=FALSE), method = "complete", members = NULL), labels=rarefied.si$Year)


simper(rarecurve_cleansam, rarefied.si$Plant, permutations = 0, trace = FALSE)
simper(rarecurve_cleansam, rarefied.si$plantyear, permutations = 0, trace = FALSE)
simper(rarecurve_cleansam, rarefied.si$Site, permutations = 0, trace = FALSE)

###WILCOXON-TEST#####
wilcox.test(shannon~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)
wilcox.test(simpson~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)
wilcox.test(richness~Plant, data = rarefied.si, exact = FALSE, correct = FALSE, conf.int = TRUE)

#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in metadata and asv table
nut_enrich_asv_table <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=3000)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(nut_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#20.4
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#12.1

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color=Type2))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 20.4%")+
  ylab("PC2- 12.1%")

ggplot(ko.coords, aes(MDS1, LogBd))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("Axis 1 values")+
  ylab("Log Bd ITS")+
  geom_smooth(method='lm')+
  stat_cor(method = "spearman", cor.coef.name="rho")


#load in asv table
asv.tbl <- read.delim("C:/Users/patty/Downloads/nut_enrich_asv_table.txt", row.names=1, header=T)
meta<-read.delim("C:/Users/patty/Downloads/nut_enrich_16S_map.txt", header=T)

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(nut_enrich_asv_table), sample=3000)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(rrarefy(t(nut_enrich_asv_table), sample = 3000))`)
t.test(larv.alph$Richness, larv.alph$Type2)


#plot richness
ggplot(larv.alph, aes(Type2, Richness, fill=Type2))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("ASV Richness")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))


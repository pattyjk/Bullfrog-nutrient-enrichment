#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in metadata and asv table
asv_table <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)

#filter out chloroplast/mitochondria
asv_table<-asv_table[-which(row.names(asv_table)=='180ed1cd4fafc390ffc300108ccf648e' | row.names(asv_table)=='7827defa226f025727dd0d1866cb5bba' | row.names(asv_table)=='94fc02cb626b227105e3b90cc5900802'),]

#look at sequencing depth
colSums(asv_table)
#2573 is the lowest good depth, lose only 3 samples

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=2573)

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
#12.2

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color=Type))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 20.4%")+
  ylab("PC2- 12.2%")

#plot Bd load as a function axis 1 values
ggplot(ko.coords, aes(MDS1, LogBd))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("Axis 1 values")+
  ylab("Log Bd ITS")+
  geom_smooth(method='lm')+
  stat_cor(method = "spearman", cor.coef.name="rho")
#rho=0.66, p=3.7e-6

###Calculate alpha diversity

#load in asv table and metadata
nut_enrich_asv_table<- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(nut_enrich_asv_table), sample=3000)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(rrarefy(t(nut_enrich_asv_table), sample = 3000))`)
t.test(larv.alph$Richness, larv.alph$Type2)
#t = 22, df = 41, p-value < 2.2e-16

larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(nut_enrich_asv_table), sample=3000), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)

t.test(larv.alph$Shannon, larv.alph$Type2)
t = 25.796, df = 41, p-value < 2.2e-16

#plot richness
ggplot(larv.alph, aes(Type, Richness, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTU Richness")+
  scale_fill_manual(values = c('#f58231', '#4363d8'))

ggplot(larv.alph, aes(Type, Shannon, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")+
  scale_fill_manual(values = c('#f58231', '#4363d8'))

#plot Bd load between pond type
ggplot(larv.alph, aes(Type, Bd_load, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  scale_y_log10()+
  ylab("Log Bd ITS Copy Number")+
  scale_fill_manual(values = c('#f58231', '#4363d8'))

#plot richness as function of Bd load
ggplot(larv.alph, aes(Richness, Bd_load, color=Type))+
  geom_point()+
  theme_bw()+
  xlab("")+
  xlab("sOTU richness")+
  ylab("Log Bd ITS Copy Number")+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  stat_cor(method = "spearman", cor.coef.name="rho")
###############################################
#analyze BioLog EcoPlate and metadata
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)
ecolog_data <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/ecolog_data_nowat.txt")

#calculate AWCD for each sample
awcd<-as.data.frame(apply(ecolog_data[,-1],2,mean))
awcd$SampleID<-row.names(awcd)
awcd<-merge(awcd, meta, by.x='SampleID', by.y='Pond', all.x=T, all.y=F)

ggplot(awcd,aes(Type,`apply(ecolog_data[, -1], 2, mean)`, fill=Type))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c('#f58231', '#4363d8'))+
  ylab("Average Well Color Development")+
  xlab("")

#calculate Shannon Diversity
library(vegan)
ecolog_data <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/ecolog_data_nowat.txt", row.names=1)
eco_shan<-as.data.frame(vegan::diversity(t(ecolog_data)))
eco_shan$Pond<-row.names(eco_shan)
names(eco_shan)<-c("Shannon", "Pond")
eco_shan$Type<-c("Reference", "Reference", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Reference", "Reference", "Reference")

ggplot(eco_shan, aes(Type, Shannon, fill=Type))+
  geom_boxplot()+
  theme_bw()+
  geom_jitter()+
  scale_fill_manual(values = c('#f58231', '#4363d8'))+
  ylab("Shanon Diversity")+
  xlab("")

t.test(eco_shan$Shannon ~ eco_shan$Type)
#no significant difference, t = -0.53257, df = 6.2581, p-value = 0.6127

scaled_eco<-scale(ecolog_data[,-1])
dist_mat <- dist(t(scaled_eco), method = 'manhattan')
hclust_avg <- hclust(dist_mat, method = 'complete')
library(ggdendro)

plot(hclust_avg)

dhc <- as.dendrogram(hclust_avg)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
 ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))+
   theme_minimal()+
   ylab("Eucliden Distance")+
   xlab("")


ggdendrogram(hclust_avg,
             rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "", y = "Euclidean distance")+
  coord_flip()

#calculate PCA
ecolog_data <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/ecolog_data_nowat.txt", row.names=1)
nut_enrich_pca<-eco_pca<-prcomp(t(ecolog_data))
pca_coords<-as.data.frame(nut_enrich_pca$x)
pca_coords$SampleID<-row.names(pca_coords)
pca_coords<-merge(pca_coords, meta, by.x='SampleID', by.y='Pond', all.y=F)

ggplot(pca_coords, aes(PC1, PC2, color=Type))+
  geom_point()+
 # geom_label()+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()

#plot carbon source development
library(reshape2)
ecolog_data <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/ecolog_data_nowat.txt")
eco_m<-melt(ecolog_data)
eco_m<-merge(eco_m, meta, by.x='variable', by.y='Pond', all.x=T, all.y=F)

#make column that removes number from carbon sources
eco_m$Carbon2<-eco_m$Carbon
eco_m$Carbon2<-gsub("[[:digit:]]", "", eco_m$Carbon2)

ggplot(eco_m, aes(Carbon2, value, fill=Type))+
  coord_flip()+
  theme_bw()+
geom_boxplot()+
  scale_fill_manual(values = c('#f58231', '#4363d8'))

#calculate Kruskal-wallis test to determine carbon sources that have different activity
library(ggpubr)
kruskal_results_eco<-as.data.frame(compare_means(value ~ Type, group.by = 'Carbon', p.adjust.method='BH', method = 'kruskal.test', data=eco_m))
#7/31 don't signifcantly differ


#####Percent inhibitory
########################

#calculate with BD database/vsearch, done in linux shell
vsearch -usearch_global rep_seqs/dna-sequences.fasta -db antiBd_db/AmphibBac_Inhibitory_2023.2r.fasta --strand plus --id 0.99 --blast6out bullfrog_out.txt

#read in results from vsearch clustering against AmphiBac database
inhibitory<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/bullfrog_out.txt", header=F)
#ASVs are V1 in df

#read in metadata and asv table
asv_table <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)

#filter out chloroplast/mitochondria
asv_table<-asv_table[-which(row.names(asv_table)=='180ed1cd4fafc390ffc300108ccf648e' | row.names(asv_table)=='7827defa226f025727dd0d1866cb5bba' | row.names(asv_table)=='94fc02cb626b227105e3b90cc5900802'),]

#rarefy table
asv_table<-rrarefy(t(asv_table), sample=2573)
asv_table<-as.data.frame(t(asv_table))

#subset ASV table to only include ASVs that match database
inhib_tb<-asv_table[row.names(asv_table) %in% inhibitory$V1,]

#calculate percentage of ASVs that are inhibitory
100*(220/3697)
#5.95%

#calculate colSums for total/inhibitory communities
total_sum<-as.data.frame(colSums(asv_table))
inhib_sum<-as.data.frame(colSums(inhib_tb))

#bind data together
inhib_tb2<-cbind(total_sum, inhib_sum)

#calculate percent inhibitory
inhib_tb2$per_inhib<-inhib_tb2$`colSums(inhib_tb)`/inhib_tb2$`colSums(asv_table)`

#add column for SampleID and merge metadata
inhib_tb2$SampleID<-row.names(inhib_tb2)
inhib_tb2<-merge(inhib_tb2, meta, by='SampleID')
inhib_tb2$per_inhib2<-inhib_tb2$per_inhib*100

#plot per species
library(ggplot2)
ggplot(inhib_tb2[-which(inhib_tb2$per_inhib2>50),], aes(Type, per_inhib2, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  scale_fill_manual(values = c('#f58231', '#4363d8'))+
  ylab("Percent Inhibitory towards Bd")

#calculate stats
TukeyHSD(aov(inhib_tb2$per_inhib ~inhib_tb2$Type))
#                       diff        lwr       upr    p adj
#Reference-Enriched 0.05785815 -0.0455522 0.1612685 0.264877


###Calculate ASVs that respond (+/-) to n-enrichment
#read in metadata and asv table
asv_table <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)

#filter out chloroplast/mitochondria
asv_table<-asv_table[-which(row.names(asv_table)=='180ed1cd4fafc390ffc300108ccf648e' | row.names(asv_table)=='7827defa226f025727dd0d1866cb5bba' | row.names(asv_table)=='94fc02cb626b227105e3b90cc5900802'),]

#rarefy table
asv_table<-rrarefy(t(asv_table), sample=2573)
asv_table<-as.data.frame(t(asv_table))

#add ASV names as a column
asv_table$OTU<-row.names(asv_table)

#reshape table for calculation
asv_m<-melt(asv_table)

#add metadata
asv_m<-merge(asv_m, meta, by.x='variable', by.y='SampleID')
asv_m<-asv_m[,c(2,3,6)]

#calculate kruskal-wallis to ID differential abundant taxa with Benjamini-hochberg corrected p-value
library(ggpubr)
kruskal_results<-as.data.frame(compare_means(value ~ Type, group.by = 'OTU', p.adjust.method='fdr', method = 'kruskal.test', data=asv_m))
length(which(kruskal_results$p.adj<0.01))
#82 differentially abundant OTUs because of N-enrichment

#read in taxonomy data for each OTU & append to Kruskal-wallis results
tax<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/taxonomy.tsv", header=T)
kruskal_results<-merge(kruskal_results, tax, by.x='OTU', by.y='Feature.ID', all.y=F)

#create a table of only significant OTUs
sig_krusk<-kruskal_results[which(kruskal_results$p.adj<0.01),]

library(plyr)
library(stringi)
library(tidyr)
#add and split taxonomy for Kruskal-Wallis test
sig_krusk<-separate(sig_krusk, Taxon, sep=';', , into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#calculate average abundance of each  for n-enrichment/reference sites
asv_means<-ddply(asv_m, c("OTU", 'Type'), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#add averages to significant kruskal wallis results
sig_krusk<-merge(sig_krusk, asv_means, by.x='OTU', by.y='OTU', all.y=F)

#change to relative abundance
sig_krusk$rel_abun<-10*(sig_krusk$mean/2573)
sig_krusk$se2<-sig_krusk$se/10

#fix genus names
write.table(sig_krusk, 'sig_krusk.txt', row.names=F, sep='\t', quote=F)
sig_krusk <- read.delim("~/GitHub/Bullfrog-nutrient-enrichment/sig_krusk.txt")

#plot mean abundance
ggplot(sig_krusk, aes(Genus, rel_abun, color=Type))+
  geom_point()+
  #facet_wrap(~Phylum)+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()+
  coord_flip()+
  ylab("Relative Abundance")+
  #facet_wrap(~Phylum)+
  xlab("")+
  geom_errorbar(aes(ymin=rel_abun-se2, ymax=rel_abun+se2, colour=Type), width=.2)

  
ggplot(sig_krusk, aes(Genus, Type, fill=rel_abun))+
  geom_tile()+
  coord_flip()+
  xlab("")+
  ylab("")+
  theme_bw()+
  scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")




###############
#picrust data
not_norm<-read.delim("Bullfrog-nutrient-enrichment/picrust_output/picrust_no_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv", row.names=1)
norm<-read.delim("Bullfrog-nutrient-enrichment/picrust_output/picrust_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv", row.names=1)

#divide table
copy_no<-not_norm/norm

#remove NAs from dividing by zero
copy_no<-na.omit(copy_no)

#calculate average
colavgs<-as.data.frame(colMeans(copy_no))

#append sample name to table
colavgs$SampleID<-row.names(colavgs)

#read in append metadata to averages
meta<-read.delim("~/GitHub/Bullfrog-nutrient-enrichment/nut_enrich_16S_map.txt", header=T)
colavgs<-merge(colavgs, meta, by='SampleID')

#plot
library(ggplot2)
ggplot(colavgs, aes(Type, `colMeans(copy_no)`, fill=Type))+
  geom_boxplot()+
  theme_bw()

#check stats
summary(aov(colavgs$`colMeans(copy_no)`~colavgs$Type))
#not significant
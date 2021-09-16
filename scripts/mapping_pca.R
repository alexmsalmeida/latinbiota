# load libraries
library(ape)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(reshape2)
library(data.table)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/uhgg-combined/")
bwa.prev = as.data.frame(read.csv("summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE))
bwa.prev = bwa.prev[,-ncol(bwa.prev)]

# load metadata
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
metadata = metadata[rownames(bwa.prev),]

# PCoA
dist.df = vegdist(bwa.prev, method="jaccard")
pcoa.data = pcoa(dist.df)
pcoa.axes = data.frame(pcoa.data$vectors)
pc1_var = round(pcoa.data$values[1,2]*100,1)
pc2_var = round(pcoa.data$values[2,2]*100,1)
bwa.df.pca = data.frame(row.names=rownames(pcoa.axes),
                        PC1=pcoa.axes[,1],
                        PC2=pcoa.axes[,2],
                        Study=metadata[rownames(pcoa.axes),"Study"],
                        State=metadata[rownames(pcoa.axes),"Health state"],
                        Continent=metadata[rownames(pcoa.axes),"Continent"],
                        Country=metadata[rownames(pcoa.axes),"Country"],
                        Age_group=metadata[rownames(pcoa.axes),"Age group"])

pca.plot = ggplot(bwa.df.pca, aes(x=PC1, y=PC2, colour=Continent)) + 
  geom_point(size=0.8, alpha=1) +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  scale_colour_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#b2450a")) + # Continents
  #scale_colour_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")) + # Countries
  ylab(paste("PC2"," (",pc2_var,"%)",sep="")) + 
  xlab(paste("PC1"," (",pc1_var,"%)",sep="")) +
  theme(legend.position="bottom", legend.box = "horizontal", legend.text=element_text(size=12),
        legend.title = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(nrow = 2, title.position="top", title.hjust = 0,
                               override.aes = list(size=2))) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))
ggsave(file="pca_euproks-above0_binary_world-adult-healthy.png", width=7, height=5, dpi=300)

# calculate % explained variance by factor
perm.res = adonis2(bwa.prev ~ Continent, metadata, method="jaccard")

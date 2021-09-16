# load libraries
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
bwa.prev = as.data.frame(read.csv("viruses-above10_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE))
bwa.prev = as.data.frame(read.csv("viruses-above0_binary_latinbiota-all.csv", row.names=1, check.names=FALSE))
bwa.prev = bwa.prev[,-ncol(bwa.prev)]
metadata = read.delim("../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
metadata = metadata[rownames(bwa.prev),]

# subset data
phage = "GUTPHAGE_37596"
variable = "Country"
samples = rownames(bwa.prev)[which(bwa.prev[,phage] > 0)]
cont.freq = data.frame(sort(table(metadata[samples,variable])/table(metadata[,variable])*100, decreasing=TRUE))

# plot barplot
bar.plot = ggplot(cont.freq, aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", alpha=0.7, colour="darkgrey", size=0.1) +
  theme_classic() +
  coord_flip() +
  ylab("% Samples") +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", 
                               "#FC8D62", "#b2450a", "#278d29", "#3b35ce")) +
  guides(fill=FALSE) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))

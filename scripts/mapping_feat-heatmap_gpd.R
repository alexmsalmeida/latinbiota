# load libraries
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(scales)
library(ggplot2)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/gpd-combined/")
bwa.data = as.data.frame(read.csv("summary_tables/viruses-above0_binary_latinbiota-all.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE))
chisq.res = read.delim("chisq_viruses_latinbiota.tsv", stringsAsFactors = FALSE)
chisq.res = chisq.res[order(chisq.res$p_value),]
taxa.all = read.delim("../../viruses/DemoVir_assignments.txt", row.names=1, stringsAsFactors = FALSE)
gpd.stats = read.delim("../../viruses/gpd-combined_stats.tsv", stringsAsFactors = FALSE)
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata = metadata[which(metadata$Study == "Latinbiota" & metadata$`Health state` == "Healthy" & metadata$`Age group` == "Adult" & metadata$Antibiotics != "Yes"),]

# subset features
feat.heat = bwa.data[metadata$Run,chisq.res[1:100,"Genome"]]
feat.heat = feat.heat[,names(sort(colSums(feat.heat), decreasing=TRUE))]

# prepare annotation
feat.taxa = data.frame(row.names=colnames(feat.heat), Family=rep(NA, ncol(feat.heat)))
feat.taxa$Family = taxa.all[match(colnames(feat.heat), rownames(taxa.all)), "Family"]
feat.taxa$Family = ifelse(is.na(feat.taxa$Family), "Unassigned", feat.taxa$Family)
taxa.colors = c("#8DD3C7", "#FFFFB3", "darkgrey")
names(taxa.colors) = sort(unique(feat.taxa$Family))
samp.geo = data.frame(Location=metadata[,"Country"], row.names=metadata[,"Run"])
geo.colors = c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")
names(geo.colors) = sort(as.vector(unique(samp.geo$Location)))
annot.colors = list(Family = taxa.colors, Location=geo.colors)

# plot heatmap
pheatmap(t(feat.heat), annotation_col = samp.geo, annotation_row = feat.taxa, annotation_colors = annot.colors, legend = FALSE,
         annotation_names_col = FALSE, annotation_names_row = FALSE,
         cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = c("grey90", "steelblue"),
         filename="heatmap_latinbiota.pdf", width=11, height=5)

# prepare barplot
bar.data = data.frame(Genome=colnames(feat.heat), Freq=colSums(feat.heat), stringsAsFactors = FALSE)
bar.data$Source = gpd.stats[match(rownames(bar.data), gpd.stats$Code),"Set"]
bar.data$Source = factor(bar.data$Source, levels=c("Latinbiota", "Shared", "GPD"))
bar.plot = ggplot(bar.data, aes(x=Genome, y=Freq, fill=Source)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  coord_flip() + 
  scale_fill_manual(values=c("brown3", "palegreen4", "steelblue"), labels=c("Latinbiota", "Shared", "Reference")) +
  scale_x_discrete(limits=rev(bar.data$Genome)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Number of samples") + 
  theme_classic() + 
  theme(legend.position="right") +
  theme(legend.title = element_text(face="bold")) +
  theme(axis.text.y = element_blank()) + 
  theme(axis.text.x = element_text(size=12), axis.ticks = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(file="barplot_latinbiota.pdf", height=5, width=3)
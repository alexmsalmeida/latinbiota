# load libraries
library(pheatmap)
library(RColorBrewer)
library(tidyr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/uhgg-combined")
bwa.data = as.data.frame(read.csv("summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE))
chisq.res = read.delim("chisq_euproks_latinbiota.tsv", stringsAsFactors = FALSE)
chisq.res = chisq.res[order(chisq.res$p_value),]
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata = metadata[which(metadata$Study == "Latinbiota" & metadata$`Health state` == "Healthy" & metadata$`Age group` == "Adult" & metadata$Antibiotics != "Yes"),]

# subset features and samples
feat.heat = bwa.data[metadata$Run,chisq.res[1:100,"Genome"]]
feat.heat = feat.heat[,names(sort(colSums(feat.heat), decreasing=TRUE))]

# prepare taxonomy
gen.stats = read.delim("../../phylogeny/uhgg-combined_stats.tsv", stringsAsFactors = FALSE, header=FALSE)
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa.proks = separate(data = gen.stats, col = V3, sep = ";", into = ranks)[,c("V1", "Order")]
colnames(taxa.proks) = c("Genome", "Order")
taxa.euks = read.delim("../../phylogeny/uhgg-euks_tax.tsv", stringsAsFactors = FALSE, header=FALSE)
colnames(taxa.euks) = c("Genome", "Order")
taxa.all = rbind(taxa.proks, taxa.euks)
rownames(taxa.all) = taxa.all$Genome

# prepare annotation
taxa.all = taxa.all[which(rownames(taxa.all) %in% colnames(feat.heat)),]
feat.taxa = data.frame(Order=taxa.all[,"Order"], row.names=rownames(taxa.all), stringsAsFactors = FALSE)
feat.taxa$Order = gsub(".*__", "", feat.taxa$Order)
keep.taxa = sort(names(sort(table(feat.taxa$Order), decreasing=TRUE)[1:10]))
feat.taxa$Order = ifelse(feat.taxa$Order %in% keep.taxa, feat.taxa$Order, "Other")
taxa.colors = c(brewer.pal(10, "Set3"), "darkgrey")
names(taxa.colors) = c(keep.taxa, "Other")
samp.geo = data.frame(Location=metadata[,"Country"], row.names=metadata[,"Run"])
geo.colors = c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")
names(geo.colors) = sort(as.vector(unique(samp.geo$Location)))
annot.colors = list(Order = taxa.colors, Location=geo.colors)

# plot heatmap
pheatmap(t(feat.heat), annotation_col = samp.geo, annotation_row = feat.taxa, annotation_colors = annot.colors, legend = FALSE,
         annotation_names_col = FALSE, annotation_names_row = FALSE,
         cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = c("grey90", "steelblue"),
         filename="heatmap_latinbiota.pdf", width=12, height=6)

# prepare barplot
bar.data = data.frame(Genome=colnames(feat.heat), Freq=colSums(feat.heat), stringsAsFactors = FALSE)
bar.data$Source = gen.stats[match(rownames(bar.data), gen.stats$V1),"V2"]
bar.plot = ggplot(bar.data, aes(x=Genome, y=Freq, fill=Source)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  coord_flip() + 
  scale_fill_manual(values=c("brown3", "palegreen4", "steelblue"), labels=c("Latinbiota", "Shared", "Reference")) +
  scale_x_discrete(limits=rev(bar.data$Genome)) +
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
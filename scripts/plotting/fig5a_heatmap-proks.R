# load libraries
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggplot2)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
bwa.data = as.data.frame(read.csv("prokaryotes/mapping/summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE))
chisq.res = read.delim("prokaryotes/mapping/chisq_euproks_latinbiota.tsv", stringsAsFactors = FALSE)
chisq.res = chisq.res[order(chisq.res$fdr, decreasing=FALSE),]
chisq.res = chisq.res[which(chisq.res$result == "Higher"),]
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata = metadata[which(metadata$Study == "Latinbiota" & metadata$`Health state` == "Healthy" & metadata$`Age group` == "Adult" & metadata$Antibiotics != "Yes"),]

# subset features and samples
feat.heat = bwa.data[metadata$Run,chisq.res[1:200,"Genome"]]
feat.heat = feat.heat[,names(sort(colSums(feat.heat), decreasing=TRUE))]

# prepare taxonomy
gen.stats = read.delim("prokaryotes/species_uhgg-latinbiota.tsv", stringsAsFactors = FALSE)
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa.proks = separate(data = gen.stats, col = Lineage, sep = ";", into = ranks)[,c("Genome", "Family")]
colnames(taxa.proks) = c("Genome", "Family")
taxa.euks = read.delim("eukaryotes/eukaryotes_tax.tsv", stringsAsFactors = FALSE, header=FALSE)
colnames(taxa.euks) = c("Genome", "Family")
taxa.all = rbind(taxa.proks, taxa.euks)
rownames(taxa.all) = taxa.all$Genome

# prepare annotation
taxa.all = taxa.all[which(rownames(taxa.all) %in% colnames(feat.heat)),]
feat.taxa = data.frame(Family=taxa.all[,"Family"], row.names=rownames(taxa.all), stringsAsFactors = FALSE)
feat.taxa$Family = gsub(".*__", "", feat.taxa$Family)
keep.taxa = sort(names(sort(table(feat.taxa$Family), decreasing=TRUE)[1:10]))
feat.taxa$Family = ifelse(feat.taxa$Family %in% keep.taxa, feat.taxa$Family, "Other")
order.taxa = c(names(sort(table(feat.taxa$Family[which(feat.taxa$Family != "Other")]), decreasing=TRUE)), "Other")
feat.taxa$Family = factor(feat.taxa$Family, levels=order.taxa)
taxa.colors = c(brewer.pal(10, "Set3"), "darkgrey")
names(taxa.colors) = order.taxa
samp.geo = data.frame(Country=metadata[,"Country"], row.names=metadata[,"Run"])
geo.colors = c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")
names(geo.colors) = sort(as.vector(unique(samp.geo$Country)))
annot.colors = list(Family = taxa.colors, Country=geo.colors)
feat.heat = feat.heat[,rownames(feat.taxa)[order(feat.taxa$Family)]]

# plot heatmap
pheatmap(t(feat.heat), annotation_col = samp.geo, annotation_row = feat.taxa, annotation_colors = annot.colors, legend = FALSE,
         annotation_names_col = FALSE, annotation_names_row = FALSE,
         cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = c("grey90", "steelblue"),
         filename="../../../Publications/2022-Latinbiota/Figures/figure5/heatmap_latinbiota.pdf", width=9, height=6)

# prepare barplot
bar.data = data.frame(Genome=colnames(feat.heat), Freq=colSums(feat.heat), stringsAsFactors = FALSE)
bar.data$Source = gen.stats[match(rownames(bar.data), gen.stats$Genome),"Source"]
bar.plot = ggplot(bar.data, aes(x=Genome, y=Freq, fill=Source)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  coord_flip() + 
  scale_fill_manual(values=c("brown3", "palegreen4", "steelblue"), 
                    labels=c("Latinbiota", "Shared", "UHGG"),
                    name="Genome source") +
  scale_x_discrete(limits=rev(bar.data$Genome)) +
  ylab("Number of samples") + 
  theme_classic() + 
  theme(legend.position="right") +
  theme(axis.text.y = element_blank()) + 
  theme(axis.text.x = element_text(size=12), axis.ticks = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(file="../../../Publications/2022-Latinbiota/Figures/figure5/barplot_latinbiota.pdf", height=7, width=3)
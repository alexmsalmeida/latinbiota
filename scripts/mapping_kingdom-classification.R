# load libraries
library(data.table)
library(ggplot2)
library(tidyverse)

# functions
read_bwa = function (x) {
  bwa.data = as.data.frame(fread(x, check.names = FALSE))
  rownames(bwa.data) = bwa.data$Genome
  bwa.data = bwa.data[,metadata$Run]
  return(bwa.data)
}

get_mapping_rates = function(x, y) {
  euproks.agg = aggregate(. ~ Kingdom, data=x, FUN=sum)
  viruses.agg = aggregate(. ~ Kingdom, data=y, FUN=sum)
  all.king = rbind(euproks.agg, viruses.agg)
  all.merged = data.frame(t(all.king[,-1]))
  colnames(all.merged) = all.king$Kingdom
  all.merged.prop = all.merged/metadata[rownames(all.merged), "Read count"]*100
  all.merged.prop$Country = metadata[rownames(all.merged), "Country"]
  all.merged.prop$Continent = metadata[rownames(all.merged), "Continent"]
  all.merged.prop = all.merged.prop[which(all.merged.prop$Country != "NA"),] 
  return(all.merged.prop)
}

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
metadata = read.delim("../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata$Country = gsub("United Republic of Tanzania", "Tanzania", metadata$Country)
rownames(metadata) = metadata$Run
bac = read.delim("../phylogeny/gtdbtk.bac120.summary.tsv", stringsAsFactors = FALSE)[,"user_genome"]
euproks = read_bwa("uhgg-combined/raw/bwa_counts-total.csv")
viruses = read_bwa("gpd-combined/raw/bwa_counts-total.csv")
viruses = viruses[,colnames(euproks)]
gpd.new = read.delim("../viruses/gpd-combined_stats.tsv", stringsAsFactors = FALSE)
gpd.new = gpd.new[which(gpd.new$Set == "Latinbiota"),1]
uhgg.new = read.delim("../phylogeny/uhgg-combined_stats.tsv", stringsAsFactors = FALSE, header=FALSE)
uhgg.new = uhgg.new[which(uhgg.new$V2 == "Latinbiota"),1]

# get mapping rates
euproks$Kingdom = ifelse(rownames(euproks) %in% bac, "Bacteria", ifelse(grepl("GC.*", rownames(euproks)), "Eukaryotes", "Archaea"))
viruses$Kingdom = "Viruses"
all.mapping = get_mapping_rates(euproks, viruses)
old.mapping = get_mapping_rates(euproks[which(!rownames(euproks) %in% uhgg.new),],
                                viruses[which(!rownames(viruses) %in% gpd.new),])

# calculate improvement
improv = (all.mapping[,c("Bacteria", "Viruses")] - old.mapping[,c("Bacteria", "Viruses")]) / old.mapping[,c("Bacteria", "Viruses")] * 100
improv$Country = all.mapping$Country
improv$Continent = all.mapping$Continent

# prepare data for plotting
all.mapping.plot = reshape2::melt(all.mapping, id.vars=c("Country", "Continent"))
all.mapping.plot = all.mapping.plot[which(all.mapping.plot$variable %in% c("Bacteria", "Viruses")),]
improv.mapping.plot = reshape2::melt(improv, id.vars=c("Country", "Continent"))
improv.mapping.plot = improv.mapping.plot[which(improv.mapping.plot$variable == "Viruses"),]

# order total mapping rate
tmp = all.mapping.plot[which(all.mapping.plot$variable == "Bacteria"),]
order.bacteria = reorder(tmp$Country, -tmp$value, FUN=median)

# plots
box.plot.all = ggplot(all.mapping.plot, aes(x=Country, y=value, fill=Continent)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~ variable, nrow=2, strip.position="right") +
  theme_bw() +
  scale_x_discrete(limits=levels(order.bacteria)) +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", 
                               "#FC8D62", "#b2450a", "#278d29", "#3b35ce")) +
  ylab("Mapped reads (%)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(strip.background=element_rect(fill="lightgrey", color=NA),
          strip.text=element_text(size=16)) +
  theme(panel.grid.minor.y = element_blank())
ggsave(file="mapping_total.png", height=6, width=15, dpi=300)

box.plot.improv = ggplot(improv.mapping.plot, aes(x=reorder(Country, -value, FUN=median), y=value, fill=Continent)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~ variable, nrow=2, strip.position="right") +
  theme_bw() +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", 
                             "#FC8D62", "#b2450a", "#278d29", "#3b35ce")) +
  coord_cartesian(ylim=c(0,100)) +
  ylab("Improvement (%)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(strip.background=element_rect(fill="lightgrey", color=NA),
        strip.text=element_text(size=16)) +
  theme(panel.grid.minor.y = element_blank())
ggsave(file="mapping_improv.png", height=6, width=15, dpi=300)

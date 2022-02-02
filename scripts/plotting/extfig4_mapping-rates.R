# load libraries
library(data.table)
library(ggplot2)
library(ggpubr)
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
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata$Country = gsub("United Republic of Tanzania", "Tanzania", metadata$Country)
metadata$Continent = ifelse(metadata$Continent == "South America" | metadata$Country == "Mexico", "Latin America", metadata$Continent)
rownames(metadata) = metadata$Run
euproks = read_bwa("prokaryotes/mapping/raw/bwa_counts-total.csv")
viruses = read_bwa("viruses/mapping/raw/bwa_counts-total.csv")
viruses = viruses[,colnames(euproks)]
gpd.new = read.delim("viruses/species_gpd-latinbiota.tsv", stringsAsFactors = FALSE)
gpd.new = gpd.new[which(gpd.new$Source == "Latinbiota"),1]
uhgg.new = read.delim("prokaryotes/species_uhgg-latinbiota.tsv", stringsAsFactors = FALSE)
bac = uhgg.new[which(grepl("Bacteria", uhgg.new$Lineage)),"Genome"]
uhgg.new = uhgg.new[which(uhgg.new$Source == "Latinbiota"),1]

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
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  ylab("Mapped reads (%)") +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=14)) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=16)) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())

box.plot.improv = ggplot(improv.mapping.plot, aes(x=reorder(Country, -value, FUN=median), y=value, fill=Continent)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~ variable, nrow=2, strip.position="right") +
  theme_bw() +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  coord_cartesian(ylim=c(0,100)) +
  ylab("Improvement (%)") +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=14)) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=16)) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())

# arrange and save
ggarrange(box.plot.all, box.plot.improv, ncol=1, nrow=2, labels = c("a", "b"), common.legend = TRUE,
          font.label = list(size = 24, color = "black", face = "bold", family = NULL), hjust=0.03, vjust=0.1)
ggsave(file="../../../Publications/2022-Latinbiota/Figures/extfig4/ExtendedData_Figure4.pdf", height=12, width=14)

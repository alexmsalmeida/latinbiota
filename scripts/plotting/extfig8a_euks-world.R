# load libraries
library(data.table)
library(matrixStats)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata$Continent = ifelse(metadata$Continent == "South America" | metadata$Country == "Mexico", "Latin America", metadata$Continent)
rownames(metadata) = metadata$Run
gen.stats = read.delim("eukaryotes/eukaryotes_tax.tsv", header=FALSE, stringsAsFactors = FALSE)
bwa.data = as.data.frame(read.csv("prokaryotes/mapping/summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE)); variable = "Continent"
euproks = bwa.data[,c(which(grepl("GC.*", colnames(bwa.data))), ncol(bwa.data))]

# group by metadata
euproks$Variable = metadata[rownames(euproks), "Continent"]
bwa.cont.group = aggregate(. ~ Variable, data=euproks, FUN=sum)
meta.counts = data.frame(table(metadata[which(metadata$Run %in% rownames(bwa.data)),variable]))

# prepare freq
euk.cont.freq = data.frame(t(bwa.cont.group[,-1]))
colnames(euk.cont.freq) = meta.counts$Var1
euk.cont.freq$Taxon = gen.stats[match(rownames(euk.cont.freq), gen.stats$V1),"V2"]
euk.cont.freq = euk.cont.freq[order(rowSums(euk.cont.freq[,-ncol(euk.cont.freq)]), decreasing=TRUE),]
euk.cont.freq$Taxon = factor(euk.cont.freq$Taxon, levels=euk.cont.freq$Taxon)
euk.cont.freq.plot = reshape2::melt(euk.cont.freq, id.vars=c("Taxon"))

# plot barplots
bar.plot.freq = ggplot(euk.cont.freq.plot, aes(x=Taxon, y=value, fill=variable)) +
  geom_bar(stat="identity", alpha=0.9, color="darkgrey", size=0.1) +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62"),
                    name="Continent") +
  ylab("Number of samples") + 
  theme_classic() +
  theme(legend.position = "top") +
  theme(strip.text = element_text(size=12, face="bold"), strip.background = element_rect(fill=NA, colour=NA)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) +
  theme(plot.margin=unit(c(t=1,r=0.5,b=1.5,l=2.5),"cm")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file="../../../Publications/2022-Latinbiota/Figures/extfig8/euks_prevalence_world-adult-healthy.pdf", width=6, height=8)

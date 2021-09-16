# load libraries
library(data.table)
library(matrixStats)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
metadata = read.delim("../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
latinbiota = metadata[which(metadata$Study == "Latinbiota" & metadata$`Age group` == "Adult" & metadata$`Health state` == "Healthy" & metadata$Antibiotics != "Yes"),"Run"]
gen.stats = read.delim("../uhgg-euks_tax.tsv", header=FALSE, stringsAsFactors = FALSE)
bwa.data = as.data.frame(read.csv("uhgg-combined/summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE)); variable = "Continent"
#bwa.data = as.data.frame(read.csv("uhgg-combined/summary_tables/euproks-above0_binary_latinbiota-all.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE)); variable = "Country"
euproks = bwa.data[,c(which(grepl("GC.*", colnames(bwa.data))), ncol(bwa.data))]

# group by metadata
bwa.cont.group = aggregate(. ~ Variable, data=euproks, FUN=sum)
meta.counts = data.frame(table(metadata[which(metadata$Run %in% rownames(bwa.data)),variable]))

# prepare freq
euk.cont.freq = data.frame(t(bwa.cont.group[,-1]))
colnames(euk.cont.freq) = meta.counts$Var1
euk.cont.freq$Taxon = gen.stats[match(rownames(euk.cont.freq), gen.stats$V1),"V2"]
euk.cont.freq = euk.cont.freq[order(rowSums(euk.cont.freq[,-ncol(euk.cont.freq)]), decreasing=TRUE),]
euk.cont.freq$Taxon = factor(euk.cont.freq$Taxon, levels=euk.cont.freq$Taxon)
euk.cont.freq.plot = reshape2::melt(euk.cont.freq, id.vars=c("Taxon"))

# prepare props
euk.cont.prop = data.frame(t(bwa.cont.group[,-1]/meta.counts$Freq*100))
colnames(euk.cont.prop) = meta.counts$Var1
euk.cont.prop$Taxon = gen.stats[match(rownames(euk.cont.prop), gen.stats$V1),"V2"]
euk.cont.prop = euk.cont.prop[order(rowMedians(as.matrix(euk.cont.prop[,-ncol(euk.cont.prop)])), decreasing=TRUE),]
euk.cont.prop$Taxon = factor(euk.cont.prop$Taxon, levels=euk.cont.prop$Taxon)
euk.cont.prop.plot = reshape2::melt(euk.cont.prop, id.vars=c("Taxon"))

# plot barplots
bar.plot.prop = ggplot(euk.cont.prop.plot, aes(x=variable, y=value, fill=variable)) +
  geom_bar(stat="identity", alpha=0.9, color="darkgrey", size=0.1) +
  facet_wrap(~ Taxon) +
  coord_flip() +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")) +
  ylab("% of samples") + 
  theme_bw() +
  guides(fill="none") +
  theme(strip.text = element_text(size=8, face="bold"), strip.background = element_rect(fill=NA, colour=NA)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.spacing.y = unit(2, "lines")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

bar.plot.freq = ggplot(euk.cont.freq.plot, aes(x=Taxon, y=value, fill=variable)) +
  geom_bar(stat="identity", alpha=0.9, color="darkgrey", size=0.1) +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#71d1ff", "#FC8D62", "#b2450a", "#3b35ce", "darkgrey")) +
  ylab("Number of samples") + 
  theme_classic() +
  guides(fill="none") +
  theme(strip.text = element_text(size=12, face="bold"), strip.background = element_rect(fill=NA, colour=NA)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) +
  theme(plot.margin=unit(c(t=1,r=0.5,b=1.5,l=2.5),"cm")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggarrange(bar.plot.freq, bar.plot.prop, widths=c(1,1.8))
ggsave(file="euks_prevalence_world-adult-healthy.png", dpi=300, width=13, height=6)
#ggsave(file="euks_prevalence_latinbiota-all.png", dpi=300, width=13, height=6)

# perform chi-squared test
n_world.euks = length(which(rowSums(euproks[,-ncol(euproks)]) > 0))
n_world.noeuks = nrow(bwa.data) - n_world.euks
n_latinbiota.euks = length(which(rowSums(euproks[latinbiota,-ncol(euproks)]) > 0))
n_latinbiota.noeuks = length(latinbiota) - n_latinbiota.euks
chisq.test(rbind(c(n_world.euks, n_world.noeuks), c(n_latinbiota.euks, n_latinbiota.noeuks)))

# investigate eukaryotic recovery vs depth
euk.freq = data.frame(Euks=rowSums(euproks[,-ncol(euproks)]))
euk.corr = merge(euk.freq, metadata, by="row.names")
corr.res = cor.test(euk.corr$Euks, euk.corr$`Read count`, method="pearson")

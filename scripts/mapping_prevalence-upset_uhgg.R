# load libraries
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
bwa.data = as.data.frame(read.csv("uhgg-combined/summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE)); variable = "Continent"
#bwa.data = as.data.frame(read.csv("uhgg-combined/summary_tables/euproks-above0_binary_latinbiota-all.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE)); variable = "Country"
metadata = read.delim("../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run

# group by cmetadata
bwa.cont.group = aggregate(. ~ Variable, data=bwa.data, FUN=sum)
meta.counts = data.frame(table(metadata[which(metadata$Run %in% rownames(bwa.data)),variable]))

# prepare upset df
upset.df = data.frame(t(bwa.cont.group[,-1]))
colnames(upset.df) = bwa.cont.group$Variable
upset.df = t(upset.df)
upset.df[upset.df/meta.counts$Freq < 0.05] = 0
upset.df[upset.df != 0] = 1
upset.df = data.frame(t(upset.df), check.names = FALSE)
upset.df = upset.df[which(rowSums(upset.df) > 0),]

# add taxonomy
taxon = "Order"
gtdb.taxa = read.delim("../phylogeny/gtdbtk.bac120.summary.tsv")
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.taxa = separate(data = gtdb.taxa, col = classification, sep = ";", into = ranks)
rownames(gtdb.taxa) = gtdb.taxa$user_genome
upset.df$Taxa = gtdb.taxa[match(rownames(upset.df), gtdb.taxa$user_genome),taxon]
upset.df$Taxa[which(is.na(upset.df$Taxa) | upset.df$Taxa == "Unassigned")] = "Unknown"
upset.df$Taxa = gsub(".*__", "", upset.df$Taxa)
keep.taxa = names(sort(table(upset.df$Taxa[which(upset.df$Taxa != "Unknown")]), decreasing=TRUE)[1:10])
upset.df$Taxa = ifelse(upset.df$Taxa %in% keep.taxa, upset.df$Taxa, "Other")
upset.df$Taxa = factor(upset.df$Taxa, levels=c(sort(keep.taxa), "Other"))
colnames(upset.df) = gsub("\\.", " ", colnames(upset.df))

# add source
gen.stats = read.delim("../phylogeny/uhgg-combined_stats.tsv", header=FALSE, stringsAsFactors = FALSE)
upset.df$Source = gen.stats[match(rownames(upset.df), gen.stats$V1),"V2"]
upset.df$Source[which(is.na(upset.df$Source))] = "Reference"
upset.df$Source = gsub("UHGG", "Reference", upset.df$Source)

# upset plot
upset(upset.df, colnames(upset.df)[-((ncol(upset.df) - 1):ncol(upset.df))], n_intersections=10, name="",
  base_annotations=list('Intersection size'=intersection_size(
    mapping=aes(fill=Source)) 
    + scale_fill_manual(values=c('Latinbiota'='brown3', 'Shared'='palegreen4', 'Reference'='steelblue'))),
  annotations = list(
    taxon=(
      ggplot(mapping=aes(fill=Taxa))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values=c(brewer.pal(10, "Set3"), "darkgrey"), name=taxon)
      + ylab('% of species'))),
  set_sizes = upset_set_size() + ylab('Number of species'),
  width_ratio=0.25,
  themes=upset_default_themes(text=element_text(size=14)))
ggsave(file="uhgg-combined/upset_euproks_world-adult-healthy.png", dpi=300, height=9, width=10)
#ggsave(file="uhgg-combined/upset_euproks_latinbiota-all.png", dpi=300, height=9, width=10)

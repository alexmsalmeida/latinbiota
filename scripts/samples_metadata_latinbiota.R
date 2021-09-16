# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# load metadata file
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
metadata = read.delim("metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata = metadata[which(metadata$Study == "Latinbiota"),]
rownames(metadata) = metadata$Run

# metadata stats
state = reshape2::melt(table(metadata$`Health state`))
state$Type = "Health state"
age = reshape2::melt(table(metadata$`Age group`))
age$Type = "Age group"
stats = rbind(state, age)
stats = stats[which(!is.na(stats$Var1)),]
stats = stats[order(stats$Type,stats$value),]
stats$Var1 = factor(stats$Var1, levels=stats$Var1)
stats$value = stats$value/nrow(metadata)*100

# plot bargraph of metadata stats
fill_colors = c(brewer.pal(5,"Blues"), brewer.pal(3,"Reds")[2:3])
meta.plot = ggplot(stats, aes(x=Type, y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.8) +
  theme_classic() +
  ylim(0,100) +
  scale_fill_manual(values = fill_colors) + 
  scale_x_discrete(limits=c("Health state", "Age group")) + 
  ylab("% Samples") + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size=12))
ggsave(file="metadata/metadata_all_barplot.pdf", height=6, width=4)

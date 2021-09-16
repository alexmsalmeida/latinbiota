# load libraries
library(UpSetR)
library(ggplot2)
library(purrr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/viral_prediction/")
#all.contigs = scan("all_contigs.txt", what="")
all.contigs = scan("filtered_contigs.txt", what="")
vs2.contigs = scan("virsorter2_contigs.txt", what="")
vibrant.contigs = scan("vibrant_contigs.txt", what="")
dvf.contigs = scan("dvf_contigs.txt", what="")
upset.df = data.frame(VirSorter2=ifelse(all.contigs %in% vs2.contigs, 1, 0),
                      VIBRANT=ifelse(all.contigs %in% vibrant.contigs, 1, 0),
                      DeepVirFinder=ifelse(all.contigs %in% dvf.contigs, 1, 0),
                      Length=as.numeric(map(strsplit(all.contigs, "_"),7)))
rownames(upset.df) = all.contigs
upset.df = upset.df[which(upset.df$Length >= 10000),]

# upset plot
pdf(file="upset_filtered.pdf", height=5, width=7, onefile=FALSE)
upset(upset.df[,!names(upset.df) == "Length"], sets = colnames(upset.df[,!names(upset.df) == "Length"]), order.by="freq", text.scale=1.5, 
      sets.bar.color = c("steelblue", "darkgreen", "darkred"))
dev.off()

# plot strength of prediction by length
upset.df$N_predictions = factor(rowSums(upset.df[,c(1:3)]))
hist.plot = ggplot(upset.df, aes(x=Length, fill=N_predictions)) + 
  geom_histogram(bins=100, colour="black", size=0.1, alpha=0.7) + 
  theme_classic() +
  xlim(5000,20000) +
  ylab("Count") + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(axis.text.x = element_text(size=10)) + 
  theme(panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("darkgreen", "steelblue", "darkred"))
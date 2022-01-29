# load libraries
library(mikropml)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# load all models
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
model.euproks.latin = read.csv("prokaryotes/mapping/machine_learning/euproks-above0_binary_latinbiota-all.csv", stringsAsFactors = FALSE)
model.euproks.latin$Database = "Prokaryotes"
model.viruses.latin = read.csv("viruses/mapping/machine_learning/viruses-above0_binary_latinbiota-all.csv", stringsAsFactors = FALSE)
model.viruses.latin$Database = "Viruses"
model.data = rbind(model.euproks.latin, model.viruses.latin)

# plot score distributions
ml.boxplots = ggplot(model.data, aes(x=method, y=Kappa, fill=Database)) +
  geom_point(alpha=0.6, position=position_jitterdodge(jitter.width=0.05, dodge.width=0.5), size=0.5, 
                    colour="darkgrey") +
  geom_boxplot(alpha=0.5, outlier.colour = NA, colour="black", width=0.5) +
  coord_flip() +
  theme_classic() +
  scale_x_discrete(limits=c("svmRadial", "rf", "glmnet"), labels=c("SVM", "Random\nForest", "Logistic\nRegression")) +
  scale_fill_manual(values=c("steelblue", "darkgreen")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  theme(legend.position="top", legend.box = "horizontal", legend.text=element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))
ggsave(file="../../../Publications/2022-Latinbiota/Figures/figure4/ml_boxplot_latinbiota.pdf", dpi=300, height=5, width=6)

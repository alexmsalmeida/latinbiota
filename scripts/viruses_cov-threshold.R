# load libraries
library(ggplot2)
library(ggrastr)
library(ggpointdensity)
library(data.table)
library(dplyr)
library(ggpubr)
library(matrixStats)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/viral_prediction/cov_threshold/")
ani99 = fread("all_vs_all_ani.tsv")
ani99$max_af = rowMaxs(as.matrix(ani99[,c("qcov", "tcov")]))
ani95 = ani99[which(ani99$pid >= 95),]

# group data
ani99.agg = ani99[ani99[, .I[which.max(max_af)], by=qname]$V1]
ani95.agg = ani95[ani95[, .I[which.max(max_af)], by=qname]$V1]

# plots
dens.plot = ggplot(ani99.agg, aes(x=max_af, y=pid)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, adjust=5) +
  ylab("ANI (%)") + 
  xlab("Coverage (%)") +
  theme_classic() +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size=16)) + 
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=14))

hist.plot = ggplot(ani95.agg, aes(x=max_af)) +
  geom_histogram() +
  ylab("Number of alignments >95% ANI") + 
  xlab("Coverage (%)") +
  geom_vline(xintercept=30, linetype="dashed") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16)) + 
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=14))

# arrange and save
ggarrange(dens.plot, hist.plot, ncol=2, nrow=1)
ggsave("cov_density.pdf", height=6, width=14)
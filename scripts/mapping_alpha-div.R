# load libraries
library(vegan)
library(stats)
library(data.table)
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
metadata = read.delim("../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
healthy.adults = metadata[which(metadata$`Age group` == "Adult" & metadata$`Health state` == "Healthy" & metadata$Antibiotics != "Yes"),"Run"]
latinbiota = metadata[which(metadata$Study == "Latinbiota" & metadata$`Age group` == "Adult" & metadata$`Health state` == "Healthy" & metadata$Antibiotics != "Yes"),"Run"]
uhgg.data = fread("uhgg-combined/summary_tables/euproks-above0_counts-transf_all.csv")
gpd.data = fread("gpd-combined/summary_tables/viruses-above0_counts-transf_all.csv")

# calculate stats
uhgg.div = diversity(uhgg.data[,-1], index="shannon")
gpd.div = diversity(gpd.data[,-1], index="shannon")
uhgg.mapped = rowSums(uhgg.data[,-1])
gpd.mapped = rowSums(gpd.data[,-1])

# prepare table
uhgg.stats = data.frame(UHGG_Diversity = uhgg.div, UHGG_Mapped = uhgg.mapped, row.names = uhgg.data$Sample)
gpd.stats = data.frame(GPD_Diversity = gpd.div, GPD_Mapped = gpd.mapped, row.names = gpd.data$Sample)
all.stats = merge(uhgg.stats, gpd.stats, by="row.names")
all.stats$Country = metadata[all.stats$Row.names, "Country"]
all.stats$Continent = metadata[all.stats$Row.names, "Continent"]
all.stats = all.stats[which(all.stats$Country != "NA"),]
all.stats$Depth = rowSums(all.stats[,c("UHGG_Mapped", "GPD_Mapped")])
all.stats.healthy = all.stats[which(all.stats$Row.names %in% healthy.adults),]
all.stats.latinbiota = all.stats[which(all.stats$Row.names %in% latinbiota),]

# visualize correlations (diversity vs depth, uhgg vs gpd)
corr.plot = ggplot(all.stats.healthy, aes(x=GPD_Diversity, y=UHGG_Diversity, colour=Continent)) +
  geom_point(size=0.2, alpha=0.8) +
  theme_classic() +
  scale_colour_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", 
                             "#FC8D62", "#b2450a")) +
  ylab("Prokaryotic diversity") +
  xlab("Viral diversity") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

# plot diversities
box.plot.all = ggplot(all.stats.healthy, aes(x=reorder(Continent, -UHGG_Diversity, FUN=median), y=UHGG_Diversity, fill=Continent)) +
  geom_point(alpha=0.6, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.5), size=0.1, 
                    colour="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  coord_flip() +
  guides(fill="none") +
  theme_bw() +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#b2450a")) +
  ylab("Alpha diversity (Shannon)") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

box.plot.latinbiota = ggplot(all.stats.latinbiota, aes(x=reorder(Country, -UHGG_Diversity, FUN=median), y=UHGG_Diversity, fill=Continent)) +
  geom_point(alpha=0.6, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.5), size=0.1, 
             colour="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  coord_flip() +
  guides(fill="none") +
  theme_bw() +
  ylab("Alpha diversity (Shannon)") +
  scale_fill_manual(values=c("#8DA0CB", "#b2450a")) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

# arrange and plot
ggarrange(corr.plot, box.plot.all, box.plot.latinbiota, nrow=1)
ggsave(file="alpha_diversity.png", dpi=300, height=5, width=14)

# statistical tests
corr.res = cor.test(all.stats.healthy$GPD_Diversity, all.stats.healthy$UHGG_Diversity, method="pearson")
wilcox.res = pairwise.wilcox.test(all.stats.healthy$UHGG_Diversity, all.stats.healthy$Country, p.adjust.method="fdr")
kruskal.res = kruskal.test(GPD_Diversity ~ Country, data = all.stats.latinbiota)

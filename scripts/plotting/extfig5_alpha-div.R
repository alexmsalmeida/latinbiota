# load libraries
library(vegan)
library(stats)
library(data.table)
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata$Continent = ifelse(metadata$Continent == "South America" | metadata$Country == "Mexico", "Latin America", metadata$Continent)
rownames(metadata) = metadata$Run
healthy.adults = metadata[which(metadata$`Age group` == "Adult" & metadata$`Health state` == "Healthy" & metadata$Antibiotics != "Yes"),"Run"]
latinbiota = metadata[which(metadata$Study == "Latinbiota" & metadata$`Age group` == "Adult" & metadata$`Health state` == "Healthy" & metadata$Antibiotics != "Yes"),"Run"]
uhgg.data = fread("prokaryotes/mapping/summary_tables/euproks-above0_counts-transf_all.csv")
gpd.data = fread("viruses/mapping/summary_tables/viruses-above0_counts-transf_all.csv")

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
all.stats$Lifestyle = metadata[all.stats$Row.names, "Lifestyle"]
all.stats = all.stats[which(all.stats$Country != "NA"),]
all.stats$Depth = rowSums(all.stats[,c("UHGG_Mapped", "GPD_Mapped")])
all.stats.healthy = all.stats[which(all.stats$Row.names %in% healthy.adults),]
all.stats.latinbiota = all.stats[which(all.stats$Row.names %in% latinbiota),]

# visualize correlations (diversity vs depth, uhgg vs gpd)
corr.plot = ggplot(all.stats.healthy, aes(x=GPD_Diversity, y=UHGG_Diversity, colour=Continent)) +
  geom_point(size=0.2, alpha=0.8) +
  theme_classic() +
  scale_colour_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
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
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  ylab("Alpha diversity") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

box.plot.latinbiota = ggplot(all.stats.latinbiota, aes(x=reorder(Country, -UHGG_Diversity, FUN=median), y=UHGG_Diversity, fill=Lifestyle)) +
  geom_point(alpha=0.6, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.5), size=0.1, 
             colour="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  coord_flip() +
  theme_bw() +
  ylab("Alpha diversity") +
  scale_fill_manual(values=c("#A6D854", "#8DA0CB")) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

# arrange and plot
ggarrange(ggarrange(corr.plot, box.plot.all, labels=c("a", "b"), legend = "bottom", common.legend = TRUE, font.label = list(size=18), vjust = 1, ncol=2, nrow=1), 
          box.plot.latinbiota, widths=c(2.3, 1), nrow=1, labels=c("", "c"), font.label = list(size=18), vjust = 1)
ggsave(file="../../../Publications/2022-Latinbiota/Figures/extfig5/ExtendedData_Figure5.pdf", height=5, width=11)

# statistical tests
corr.res = cor.test(all.stats.healthy$GPD_Diversity, all.stats.healthy$UHGG_Diversity, method="pearson")
wilcox.res = pairwise.wilcox.test(all.stats.healthy$UHGG_Diversity, all.stats.healthy$Country, p.adjust.method="fdr")
kruskal.res = kruskal.test(GPD_Diversity ~ Country, data = all.stats.latinbiota)

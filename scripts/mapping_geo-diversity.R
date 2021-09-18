# load libraries
library(vegan)
library(ggplot2)
library(cepiigeodist)
library(countrycode)
library(metagMisc)

# function
scale_data = function(x) {(x-min(x))/(max(x)-min(x))}

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/uhgg-combined/")
bwa.prev = as.data.frame(read.csv("summary_tables/euproks-above0_binary_world-adult-healthy.csv", row.names=1, check.names=FALSE))
bwa.prev = bwa.prev[,-ncol(bwa.prev)]

# load metadata
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
metadata = metadata[rownames(bwa.prev),]

# pairwise distances
dist.df = vegdist(bwa.prev, method="jaccard")
dist.pairwise = data.frame(dist2list(dist.df, tri = TRUE), stringsAsFactors = FALSE)
colnames(dist.pairwise) = c("Sample1", "Sample2", "Distance")

# add countries
dist.pairwise$Country1 = metadata[match(dist.pairwise$Sample1, metadata$Run),"Country"]
dist.pairwise$Country2 = metadata[match(dist.pairwise$Sample2, metadata$Run),"Country"]
dist.pairwise$Continent1 = metadata[match(dist.pairwise$Sample1, metadata$Run),"Continent"]
dist.pairwise$Continent2 = metadata[match(dist.pairwise$Sample2, metadata$Run),"Continent"]
dist.pairwise$Country_res = ifelse(dist.pairwise$Country1 == dist.pairwise$Country2, "Same country", "Different country")
dist.pairwise$Continent_res = ifelse(dist.pairwise$Continent1 == dist.pairwise$Continent2, "Same continent", "Different continent")

# add geo dist
geo.dists = as.data.frame(dist_cepii)
geo.dists$Codes = paste(geo.dists$iso_o, "_", geo.dists$iso_d, sep = "")
dist.geo = dist.pairwise[which(dist.pairwise$Country_res == "Different country"),]
dist.geo$Code1 = countrycode(sourcevar = dist.geo[,"Country1"], origin = "country.name", destination = "iso3c")
dist.geo$Code2 = countrycode(sourcevar = dist.geo[,"Country2"], origin = "country.name", destination = "iso3c")
dist.geo$Codes = paste(dist.geo$Code1, "_", dist.geo$Code2, sep = "")
dist.geo = merge(dist.geo, geo.dists, by="Codes")
dist.geo$Dist_scaled = scale_data(dist.geo$dist)

# plot data
plot.df = dist.pairwise[which(dist.pairwise$Continent_res == "Same continent"),]
box.plot = ggplot(plot.df, aes(x=Continent1, y=Distance, fill=Country_res)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  theme_bw() +
  ylab("Beta diversity (Jaccard distance)") +
  scale_x_discrete(limits=c("North America", "Asia", "South America", "Africa", "Europe")) +
  scale_fill_manual(values=c("steelblue", "tomato")) +
  facet_wrap(~ Country_res, nrow=2, strip.position = "right") +
  coord_flip() +
  guides(fill="none") +
  theme(strip.background = element_blank(), strip.text = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12)) +
  theme(panel.grid.minor = element_blank())
ggsave(box.plot, file="intra_vs_international_div.png", dpi=300, height=5, width=7)

# statistical tests
test.diff = dist.pairwise[which(dist.pairwise$Continent_res == "Same continent" & dist.pairwise$Country_res == "Different country"),]
kruskal.diff = kruskal.test(Distance ~ Continent1, data = test.diff)
wilcox.diff = pairwise.wilcox.test(test.diff$Distance, test.diff$Continent1, p.adjust.method="fdr")

test.same = dist.pairwise[which(dist.pairwise$Continent_res == "Same continent" & dist.pairwise$Country_res == "Same country"),]
kruskal.same = kruskal.test(Distance ~ Continent1, data = test.same)
wilcox.same = pairwise.wilcox.test(test.same$Distance, test.same$Continent1, p.adjust.method="fdr")

# plot correlation
corr.data = dist.geo[sample(rownames(dist.geo), 10000),]
corr.plot = ggplot(corr.data, aes(x=Dist_scaled, y=Distance)) +
  geom_point(alpha=0.5, size=0.01) +
  geom_smooth(method="lm") +
  theme_classic() +
  ylab("Beta diversity") +
  xlab("Geographic distance") +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12)) +
  theme(panel.grid.minor = element_blank())
ggsave(corr.plot, file="geo_vs_beta_dist.png", dpi=300, height=3, width=11)
corr.res = cor.test(dist.geo$Distance, dist.geo$Dist_scaled, method="pearson")

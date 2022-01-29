# load libraries
library(maptools)
library(ggplot2)
library(grid)
library(rgdal)
library(rgeos)
library(gpclib)
library(rgdal)
library(maptools)
library(mapproj)

# load input
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
metadata = metadata[which(metadata$Study == "Latinbiota"),]
countries = unique(metadata$Country)
countries = countries[countries != "NA"]
ddf = data.frame(matrix(NA, ncol=2, nrow=length(countries)), row.names=countries)
colnames(ddf) = c("country", "samples")
ddf$country = rownames(ddf)

# count metagenomes per country
for (c in ddf$country){
  samples = rownames(metadata)[which(metadata$Country == c)]
  ddf[c,2] = length(samples)
}
ddf$samples_class = NA
ddf$samples_class[which(ddf$samples >=25)] = "25-50"
ddf$samples_class[which(ddf$samples >50)] = "51-100"
ddf$samples_class[which(ddf$samples > 100)] = ">100"
ddf$samples_class[which(is.na(ddf$samples_class))] = "<25"
ddf$samples_class = factor(ddf$samples_class, levels=c("<25", "25-50", "51-100", ">100"))

# edit world map template
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
wrld = fortify(wrld_simpl, region="id")
wrld = subset(wrld, id != "Antarctica")

# plot map
map.plot = ggplot() + 
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="white", color="darkgrey", size=0.1) + 
  geom_map(data=ddf, map=wrld, aes(map_id=country, fill=samples_class),  color="black", size=0.2, alpha=0.7) + 
  scale_fill_manual(values=c("darkolivegreen1", "chartreuse3", "darkolivegreen4", "darkgreen"), name="Number of samples") +
  coord_map() + 
  labs(x="", y="") +
  coord_cartesian(ylim=c(-55,30), xlim=c(-115,-35)) +
  theme(panel.background = element_rect(fill = "azure"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=0.2, fill=NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")
ggsave("../../../Publications/2022-Latinbiota/Figures/figure1/map_latinbiota.pdf", height=6, width=5)

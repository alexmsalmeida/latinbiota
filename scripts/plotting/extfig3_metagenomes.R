# load libraries
library(RColorBrewer)
library(ggpubr)
library(maptools)
library(data.table)
library(ggplot2)
library(reshape2)
library(splitstackshape)
library(RColorBrewer)
library(grid)
library(rgdal)
library(rgeos)
library(gpclib)
library(rgdal)
library(maptools)
library(mapproj)
library(countrycode)

# load input
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata$Continent = ifelse(metadata$Continent == "South America" | metadata$Country == "Mexico", "Latin America", metadata$Continent)

# prepare dataset
ddf = data.frame(table(metadata$Country))
colnames(ddf) = c("country", "samples")
ddf$samples_class = NA
ddf$samples_class[which(ddf$samples >=10)] = "10-99"
ddf$samples_class[which(ddf$samples >=100)] = "100-1000"
ddf$samples_class[which(ddf$samples > 1000)] = ">1000"
ddf$samples_class[which(is.na(ddf$samples_class))] = "<10"
ddf$samples_class = factor(ddf$samples_class, levels=c("<10", "10-99", "100-1000", ">1000"))

# edit world map template
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
wrld = fortify(wrld_simpl, region="id")
wrld = subset(wrld, id != "Antarctica")

# plot map
wrld.plot = ggplot() + 
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="lightgray", color="white", size=0.1) + 
  geom_map(data=ddf, map=wrld, aes(map_id=country, fill=samples_class),  color="white", size=0.2, alpha=0.7) + 
  scale_fill_manual(values=c("steelblue1", "steelblue", "darkblue", "darkgreen"), name="Number of metagenomes") + 
  coord_map() + 
  labs(x="", y="") + 
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(colour = "black", size=0.5, fill=NA), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")

# plot continent proportion
ddf.cont = data.frame(table(metadata$Continent))
colnames(ddf.cont) = c("continent", "samples")
ddf.cont = ddf.cont[which(ddf.cont$continent != "NA"),]
ddf.cont$proportion = ddf.cont$samples / sum(ddf.cont$samples) * 100
cont.prop = ggplot(ddf.cont, aes(x=continent, y=proportion, fill=continent)) +
  geom_bar(stat="identity", colour="darkgrey", alpha=0.8) +
  theme_classic() +
  ylab("% Metagenomes") +
  scale_x_discrete(limits=ddf.cont[order(ddf.cont$samples, decreasing=TRUE),"continent"],
                   label = gsub(" ", "\n", ddf.cont[order(ddf.cont$samples, decreasing=TRUE),"continent"])) +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  guides(fill="none") +
  theme(axis.text.x = element_text(size=9)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_blank())

# arrange plot with map
ggarrange(wrld.plot, cont.prop, nrow=1, widths=c(1.8,1), labels=c("a", "b"), font.label = list(size = 18, color = "black", face = "bold", family = NULL))
ggsave(file="../../../Publications/2022-Latinbiota/Figures/extfig3/ExtendedData_Figure3.pdf", height=5, width=10)
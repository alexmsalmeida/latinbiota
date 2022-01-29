# load libraries
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(ggpubr)

# load input files
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
amr.stats = read.delim("prokaryotes/amr_profiles/amr_classes.tsv", check.names=FALSE, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(amr.stats) = c("Sample", "Class", "Counts")
amr.counts = read.delim("prokaryotes/amr_profiles/amr_counts.tsv", check.names=FALSE, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(amr.counts) = c("Sample", "Gene", "Counts", "Gene_length")

# load metadata
metadata.all = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
metadata.all = metadata.all[which(metadata.all$`Health state` == "Healthy" & metadata.all$`Age group` == "Adult" & metadata.all$Antibiotics != "Yes"),]
metadata.all$Continent = ifelse(metadata.all$Continent == "South America" | metadata.all$Country == "Mexico", "Latin America", metadata.all$Continent)
metadata.all$Country = gsub("United Republic of Tanzania", "Tanzania", metadata.all$Country)
rownames(metadata.all) = metadata.all$Run
countries.keep = names(which(table(metadata.all$Country) >= 10))
metadata = metadata.all[unique(amr.stats$Sample),]
countries = unique(metadata$Country)
countries = countries[countries != "NA"]

# calculate proportion of AMR reads per metagenome
amr.counts.agg = aggregate(Counts ~ Sample, data=amr.counts, FUN=sum)
rownames(amr.counts.agg) = amr.counts.agg$Sample
amr.counts.agg = merge(amr.counts.agg, metadata, by="row.names")
amr.counts.agg$AMR_prop = amr.counts.agg$Counts/amr.counts.agg$`Read count`*100

# calculate proportion of samples with AMR
amr.samples = data.frame(table(amr.counts.agg$Country))
total.samples = data.frame(table(metadata.all$Country))
amr.samples.prop = merge(amr.samples, total.samples, by="Var1")
amr.samples.prop$Proportion = amr.samples.prop$Freq.x/amr.samples.prop$Freq.y*100
amr.samples.prop$Continent = metadata[match(amr.samples.prop$Var1, metadata$Country),"Continent"]
order.country = as.vector(amr.samples.prop[order(amr.samples.prop$Proportion, decreasing=TRUE),"Var1"])
order.country = order.country[which(order.country %in% countries.keep)]
amr.samples.prop.cont = aggregate(list(amr.samples.prop$Freq.x, amr.samples.prop$Freq.y), by = list(amr.samples.prop$Continent), FUN=sum)
colnames(amr.samples.prop.cont) = c("Continent", "AMR", "Total")
amr.samples.prop.cont$Prop = amr.samples.prop.cont$AMR/amr.samples.prop.cont$Total*100
amr.samples.prop.cont = amr.samples.prop.cont[which(amr.samples.prop.cont$Total > 100),]

# find most prevalent AMR classes
classes.counts = aggregate(Counts ~ Class, data=amr.stats, FUN=sum)
classes.counts.notother = classes.counts[which(classes.counts$Class != "Other"),]
classes.keep = classes.counts.notother[order(classes.counts.notother$Counts, decreasing=TRUE)[1:10],"Class"]
amr.stats$Class = ifelse(amr.stats$Class %in% classes.keep, amr.stats$Class, "Other")
amr.stats$Class = factor(amr.stats$Class, levels=c(sort(classes.keep), "Other"))
amr.stats$Country = metadata[match(amr.stats$Sample, metadata$Run), "Country"]
amr.stats.agg = aggregate(Counts ~ Class+Country, data=amr.stats, FUN=sum)
amr.stats.agg$Continent = metadata[match(amr.stats.agg$Country, metadata$Country), "Continent"]

# plot proportion of samples with AMR per continent
amr.prop.cont = ggplot(amr.samples.prop.cont, aes(x=reorder(Continent, -Prop), y=Prop, fill=Continent)) +
  geom_bar(stat="identity", alpha=0.7) +
  theme_classic() +
  scale_x_discrete(labels=c("Latin\nAmerica", "Africa", "Asia", "North\nAmerica", "Europe")) +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  ylab("% Samples with AMR genes") +
  ylim(0,100) +
  guides(fill = guide_legend(title.position = "top")) +
  theme(legend.position="bottom") +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())

# plot proportion of samples with AMR per country
amr.prop = ggplot(amr.samples.prop, aes(x=Var1, y=Proportion, fill=Continent)) +
  geom_bar(stat="identity", alpha=0.7) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#b2450a", "#8DA0CB", "#FC8D62")) +
  scale_x_discrete(limits=order.country) +
  ylab("% Samples with AMR genes") +
  xlab("\n") +
  guides(fill = guide_legend(title.position = "top")) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())
amr.prop.legend = as_ggplot(get_legend(amr.prop))

# plot classes per country
stack.cont = ggplot(amr.stats.agg, aes(x=Country, y=Counts, fill=Class)) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=scales::percent_format(suffix="")) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#D9D9D9",
                             "#FDB462","#B3DE69","#FCCDE5","#80B1D3","#BC80BD","darkgrey")) +
  scale_x_discrete(limits=order.country) +
  ylab("% Detected AMR genes") +
  guides(fill = guide_legend(title.position = "top")) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(panel.grid.minor.y = element_blank())
amr.classes.legend = as_ggplot(get_legend(stack.cont))

# arrange and save plot
ggarrange(ggarrange(amr.prop.cont, amr.prop, stack.cont, nrow=1, legend = "none", 
                    widths=c(0.9,1.5,1), labels = c("a", "b", "c"), hjust = c(0,0,1), font.label = list(size=20)), 
          ggarrange(amr.prop.legend, amr.classes.legend, nrow=1, widths=c(1,1.5)), nrow=2, heights=c(5,1))
ggsave("../../../Publications/2022-Latinbiota/Figures/figure6/Figure6.pdf", height=7, width=13)

# load libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/viral_prediction/")
demovir = read.delim("DemoVir_assignments.txt")
stats = read.delim("gpd-combined_stats.tsv", stringsAsFactors = FALSE)
novel.genomes = stats$Code[which(stats$Set == "Latinbiota")]
shared.genomes = stats$Code[which(stats$Set == "Shared")]
demovir$Set = ifelse(demovir$Sequence_ID %in% novel.genomes, "Latinbiota", ifelse(demovir$Sequence_ID %in% shared.genomes, "Shared", "GPD"))
demovir.trf = as.data.frame(table(demovir[,c("Order", "Family", "Set")]))
demovir.trf = demovir.trf[which(demovir.trf$Freq > 10 & demovir.trf$Order == "Caudovirales"),]
demovir.trf$Family = factor(demovir.trf$Family, levels=c("Siphoviridae", "Myoviridae", "Podoviridae", "Unassigned"))

# calculate proportions
demovir.prop = demovir.trf %>%
  group_by(Family) %>%
  mutate(Total= sum(Freq)) %>%
  group_by(Set, .add=TRUE) %>%
  mutate(Prop=Freq/Total)
demovir.prop$Set = factor(demovir.prop$Set, levels=c("Latinbiota", "Shared", "GPD"))

# plot taxonomy
tax.plot = ggplot(demovir.prop, aes(x=Family, fill=Set, y=Freq)) +
  geom_bar(stat="identity", colour="darkgrey", alpha=0.8) +
  geom_text(aes(label=ifelse(Set=="Latinbiota", paste0(round(Prop*100,1),"%"), ""), hjust="left"), 
            position = position_stack(), hjust=-0.1) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c("darkolivegreen2", "tomato", "steelblue")) +
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.1))) +
  ylab("Number of viral clusters") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        legend.position="top")
ggsave("gpd-combined_tax.pdf", height=6, width=8)
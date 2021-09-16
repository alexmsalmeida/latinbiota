# load libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(RColorBrewer)

# load files
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mashdiff/")
status = read.delim("../../Project_resource/species_status.tab", check.names = FALSE)
prev = read.delim("uhgg_frequency.tsv", header=FALSE)
colnames(prev) = c("Genome", "Counts")
clst.data = merge(status, prev, by="Genome")

# parse taxonomy
gtdb.taxa = read.delim("../../Project_resource/species_taxonomy.tab")
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.taxa = separate(data = gtdb.taxa, col = GTDB_lineage, sep = ";", into = ranks)
rownames(gtdb.taxa) = gtdb.taxa$Genome
clst.nr = merge(clst.data, gtdb.taxa, by="Genome")
clst.nr$Taxon = NA
for (n in 1:nrow(clst.nr)){
  for (r in ranks){
    if(!grepl("__$", clst.nr[n,r])){
      clst.nr[n,"Taxon"] = clst.nr[n,r]
    }
  }
}

# barplot by type
order.nr = clst.nr$Genome[order(clst.nr[,"Counts"], decreasing=TRUE)]
values.nr = gsub("_", " ", gsub(".*__","",clst.nr$Taxon[order(clst.nr[,"Counts"], decreasing=TRUE)]))
bar.plot = ggplot(clst.nr, aes(x=Genome, y=Counts, fill=Status)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  coord_flip() + 
  scale_x_discrete(limits=order.nr[1:25], labels=values.nr[1:25]) + 
  scale_fill_manual(values=c("steelblue", "darkred", "palegreen4"), name = "Species status",
                    labels=c("Cultured (gut)", "Cultured (other)", "Uncultured")) + 
  ylab("Number of genomes") + 
  theme_classic() + 
  theme(legend.position="top") +
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave("uhgg_frequency.tiff", height=6, width=8, dpi=300)

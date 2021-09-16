# load libraries
library(phytools)
library(ggplot2)
library(ggpubr)

# load files
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/phylogeny/")
tre = read.tree("gtdbtk.bac120.user_msa.fasta.treefile") # load tree
tax = read.delim("itol_gtdb-layer.txt", skip=4L, sep=",", header=FALSE, stringsAsFactors = FALSE)
gset = read.delim("itol_genome-set.txt", skip=4L, sep=",", header=FALSE, stringsAsFactors = FALSE)
df = merge(tax, gset, by="V1")
colnames(df) = c("Genome", "Colour", "Phylum", "iTOL", "Colour2", "Status")
tax.known = df[which(df$Status != "Latinbiota"),]
tax.novel = df[which(df$Status == "Latinbiota"),]

# prepare file
all.phyla = c(as.vector(tax.known$Phylum), as.vector(tax.novel$Phylum))
phyla = unique(all.phyla[which(!is.na(all.phyla))])
pd = data.frame(matrix(0, ncol=3, nrow=length(phyla)+1))
colnames(pd) = c("Phylum", "Known_PD", "Known+Novel_PD")

# calculate phylogenetic diversity (PD)
for (n in 1:length(phyla)){
  phylum = as.character(phyla[n])
  species.hgr = as.vector(tax.known[which(tax.known$Phylum == phylum),"Genome"])
  if (length(species.hgr) > 0){
    subs.tre.hgr = drop.tip(tre, tre$tip.label[-match(species.hgr, tre$tip.label)])
    pd[n,2] = sum(subs.tre.hgr$edge.length)
  } 
  species.all = as.vector(df[which(df$Phylum == phylum),"Genome"])
  subs.tre.all = drop.tip(tre, tre$tip.label[-match(species.all, tre$tip.label)])
  pd[n,1] = phylum
  pd[n,3] = sum(subs.tre.all$edge.length)
}

# calculate improvement
tre.hgr = drop.tip(tre, tre$tip.label[-match(tax.known$Genome, tre$tip.label)])
pd[length(phyla)+1,] = c("Total", sum(tre.hgr$edge.length), sum(tre$edge.length))
pd$Proportion = ((as.numeric(pd$`Known+Novel_PD`)-as.numeric(pd$Known_PD))
                 /as.numeric(pd$`Known+Novel_PD`)*100) # proportion
pd$Improvement = (as.numeric(pd$`Known+Novel_PD`)-as.numeric(pd$Known_PD)) # difference
pd = pd[which(pd$Improvement > 0),]

# indicate total improvement
total.improv = pd[nrow(pd), "Improvement"]/as.numeric(pd[nrow(pd), "Known_PD"])*100
cat("Novel species provide a", total.improv, "% increase in phylogenetic diversity")

# reorder phylum
pd.fi = pd[order(pd$Proportion[-nrow(pd)], decreasing=TRUE),]
pd.fi$Phylum = factor(pd.fi$Phylum, levels = pd.fi$Phylum)

# define colors
color.phyla=df$Colour
names(color.phyla) = df$Phylum

# plot phylogenetic diversity
prop.total = ggplot(pd.fi, aes(x=Phylum, y=Proportion, fill=Phylum)) +
  geom_bar(stat="identity", alpha=0.75, size=0.2) +
  theme_classic() +
  scale_fill_manual(values=color.phyla) +
  ylab("Proportion of the diversity provided\nby novel species (%)") +
  guides(fill=FALSE) +
  coord_flip() + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12))

prop.improv = ggplot(pd.fi, aes(x=Phylum, y=Improvement, fill=Phylum)) +
  geom_bar(stat="identity", alpha=0.75, size=0.2) +
  theme_classic() +
  scale_fill_manual(values=color.phyla) +
  ylab("Increase in diversity provided by the\nnovel species (total branch length)") +
  guides(fill=FALSE) +
  coord_flip() + 
  #scale_y_reverse() 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12))

ggarrange(prop.total, prop.improv, ncol=1, nrow=2)
ggsave(filename="phylo_diversity.tiff", dpi=300, height=8, width=5)

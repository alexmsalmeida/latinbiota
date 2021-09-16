# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mashdiff/")
dnadiff.uhgg = read.delim("dnadiff_vs_uhgg.tsv", header=FALSE, stringsAsFactors = FALSE)
colnames(dnadiff.uhgg) = c("query", "reference", "ref_len", "ref_aln", "qry_len", "qry_aln", "ani", "mash")
dnadiff.uhgg$reference = gsub(".fa", "", dnadiff.uhgg$reference)
checkm.data = read.csv("all_checkm.csv")
checkm.data$genome = gsub(".fa", "", checkm.data$genome)
genome2species = read.delim("../../Project_resource/genomes-all_metadata.tsv", stringsAsFactors = FALSE)[,c("Genome","Species_rep")]

# prepare dataset
dnadiff.uhgg$UHGG_species = genome2species$Species_rep[match(dnadiff.uhgg$reference, genome2species$Genome)]
dnadiff.uhgg[,c("UHGG_comp", "UHGG_cont")] = checkm.data[match(dnadiff.uhgg$UHGG_species, checkm.data$genome), 
                                                         c("completeness", "contamination")]
dnadiff.uhgg$UHGG_qs = dnadiff.uhgg$UHGG_comp - 5*dnadiff.uhgg$UHGG_cont
dnadiff.uhgg[,c("MAGs_comp", "MAGs_cont")] = checkm.data[match(dnadiff.uhgg$query, checkm.data$genome), 
                                                         c("completeness", "contamination")]
dnadiff.uhgg$MAGs_qs = dnadiff.uhgg$MAGs_comp - 5*dnadiff.uhgg$MAGs_cont
dnadiff.uhgg$Match = ifelse(max(dnadiff.uhgg$ref_aln, dnadiff.uhgg$qry_aln) >= 80 & dnadiff.uhgg$ani >= 99, "Strain match",
                            ifelse(max(dnadiff.uhgg$ref_aln, dnadiff.uhgg$qry_aln) >= 30 & dnadiff.uhgg$ani >= 95, "Species match",
                                   "No match"))

# dnadiff scatterplot
ani_mash = dnadiff.uhgg[which(dnadiff.uhgg$mash < 0.1),]
ani_mash$mash = 100-ani_mash$mash*100
ani_mash$min_af = apply(ani_mash[,c("ref_aln", "qry_aln")], 1, min)
ani_mash$af_class = ifelse(ani_mash$min_af >= 60, ">=60%", ifelse(ani_mash$min_af >= 30, "30-60%", "<30%"))
dnadiff.plot = ggplot(ani_mash, aes(x=mash, y=ani, colour=af_class)) +
  geom_point(size=0.5, alpha=0.4) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_colour_manual(values=c("steelblue", "darkred"), name="Aligned fraction") +
  theme_classic() +
  xlim(90,100) +
  ylim(90,100) +
  ylab("ANImf") +
  xlab("Mash") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))
ggsave("animf_vs_mash.pdf", height=5, width=6)

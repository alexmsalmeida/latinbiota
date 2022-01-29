# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
dnadiff.uhgg = read.delim("prokaryotes/mags_vs_uhgg-all_besthits_filt.tsv", header=TRUE, stringsAsFactors = FALSE)
genome2species = read.delim("../Project_resource/genomes-all_metadata.tsv", stringsAsFactors = FALSE)[,c("Genome","Species_rep")]
gtdb.tk = read.delim("../Project_resource/species_taxonomy.tab")
colnames(gtdb.tk)[1] = "UHGG_species"

# prepare dataset
dnadiff.uhgg$UHGG_species = genome2species$Species_rep[match(dnadiff.uhgg$ref, genome2species$Genome)]
dnadiff.uhgg = merge(dnadiff.uhgg, gtdb.tk, by="UHGG_species")
dnadiff.uhgg$af_max = apply(dnadiff.uhgg[,c("ref_cov", "qry_cov")], 1, max)
dnadiff.uhgg$Match = ifelse(dnadiff.uhgg$af_max >= 80 & dnadiff.uhgg$ani >= 99, "Strain match",
                            ifelse(dnadiff.uhgg$af_max >= 30 & dnadiff.uhgg$ani >= 95, "Species match",
                                   "No match"))

# dnadiff scatterplot
dnadiff.plot = ggplot(dnadiff.uhgg, aes(x=qry_cov, y=ref_cov, colour=Match)) +
  geom_point(size=0.5, alpha=0.4) +
  theme_classic() +
  scale_colour_manual(values=c("tomato", "steelblue", "olivedrab3"), name="Match level") +
  ylab("% Reference aligned (UHGG)") +
  xlab("% Query aligned (Latinbiota)") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

dens.plot = ggplot(dnadiff.uhgg, aes(x=qry_cov, fill=Match)) +
  geom_density(alpha=0.8, size=0.1) +
  theme_classic() +
  scale_fill_manual(values=c("tomato", "steelblue", "olivedrab3"), name="Match level") +
  ylab("Density") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=12))

bar.df = data.frame(table(dnadiff.uhgg$Match))
bar.df$Match = "MAGs"
bar.df$Var1 = factor(bar.df$Var1, levels=c("Strain match", "Species match", "No match"))
bar.plot = ggplot(bar.df, aes(x=Match, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.8, size=0.1) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values=rev(c("tomato", "steelblue", "olivedrab3")), name="Match level") +
  scale_y_continuous(breaks=c(0,4000,8000,12000,16000)) +
  ylab("Number of MAGs") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# arrange plot and save
ggarrange(ggarrange(dens.plot, dnadiff.plot, ncol=1, nrow=2, heights=c(1,4), align="v", legend = FALSE), 
          bar.plot, ncol=1, nrow=2, heights=c(4,1), common.legend = TRUE)
ggsave("../../../Publications/2022-Latinbiota/Figures/extfig2/dnadiff_vs_uhgg.pdf", height=8, width=6)

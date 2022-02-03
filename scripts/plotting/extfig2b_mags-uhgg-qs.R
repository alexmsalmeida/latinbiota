# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
dnadiff.uhgg = read.delim("prokaryotes/mags_vs_uhgg.tsv", header=TRUE, stringsAsFactors = FALSE)
checkm.stats = read.delim("prokaryotes/genomes_uhgg-latinbiota.tsv", stringsAsFactors = FALSE)
genome2species = read.delim("http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv", stringsAsFactors = FALSE)[,c("Genome", "Genome_type", "Species_rep")]

# prepare dataset
dnadiff.uhgg$UHGG_species = genome2species$Species_rep[match(dnadiff.uhgg$ref, genome2species$Genome)]
dnadiff.uhgg[,c("UHGG_comp", "UHGG_cont")] = checkm.stats[match(dnadiff.uhgg$UHGG_species, checkm.stats$Genome), 
                                                         c("Completeness", "Contamination")]
dnadiff.uhgg$UHGG_qs = dnadiff.uhgg$UHGG_comp - 5*dnadiff.uhgg$UHGG_cont
dnadiff.uhgg[,c("MAGs_comp", "MAGs_cont")] = checkm.stats[match(dnadiff.uhgg$qry, checkm.stats$Genome), 
                                                         c("Completeness", "Contamination")]
dnadiff.uhgg$MAGs_qs = dnadiff.uhgg$MAGs_comp - 5*dnadiff.uhgg$MAGs_cont
dnadiff.uhgg$QS_result = ifelse(dnadiff.uhgg$MAGs_qs > dnadiff.uhgg$UHGG_qs, "Latinbiota better", "UHGG better")
dnadiff.uhgg$QS_ratio = dnadiff.uhgg$MAGs_qs/dnadiff.uhgg$UHGG_qs
dnadiff.uhgg$af_max = apply(dnadiff.uhgg[,c("ref_cov", "qry_cov")], 1, max)
dnadiff.uhgg = dnadiff.uhgg[which(dnadiff.uhgg$af_max >= 30 & dnadiff.uhgg$ani >= 95),]

# subset data
qs.df = aggregate(QS_ratio ~ UHGG_species, data=dnadiff.uhgg, FUN=max)
for (row in 1:nrow(qs.df)){
  qs.df[row,c("qry", "MAGs_comp", "MAGs_qs", "UHGG_comp", "UHGG_qs")] = dnadiff.uhgg[which(dnadiff.uhgg$UHGG_species == qs.df[row,"UHGG_species"] &
                                                       dnadiff.uhgg$QS_ratio == qs.df[row, "QS_ratio"]),
                                                  c("qry", "MAGs_comp", "MAGs_qs", "UHGG_comp", "UHGG_qs")]
}
qs.df$QS_result = ifelse(qs.df$MAGs_qs > qs.df$UHGG_qs, "Latinbiota better", ifelse(qs.df$MAGs_qs < qs.df$UHGG_qs, "UHGG better", "Equal"))

# QS scatterplot
qs.plot = ggplot(qs.df, aes(x=MAGs_qs, y=UHGG_qs, colour=QS_result)) +
  geom_point(size=0.8, alpha=0.6) +
  theme_classic() +
  scale_colour_manual(values=c("grey", "tomato", "steelblue"), name="QS comparison") +
  ylab("UHGG Quality Score") +
  xlab("Latinbiota MAGs Quality Score") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

bar.qs = data.frame(table(qs.df$QS_result))
bar.qs$QS = "Species matches"
bar.qs$Var1 = factor(bar.qs$Var1, levels=c("UHGG better", "Latinbiota better", "Equal"))
bar.qs.plot = ggplot(bar.qs, aes(x=QS, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.8, size=0.1) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values=rev(c("grey", "tomato", "steelblue")), name="QS comparison") +
  scale_y_continuous(breaks=c(0,400,800,1200)) +
  ylab("Number of species") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# arrange and save plot
ggarrange(qs.plot, bar.qs.plot, ncol=1, nrow=2, heights=c(4,1), common.legend = TRUE)
ggsave(file="../../../Publications/2022-Latinbiota/Figures/extfig2/qs_comparison.pdf", height=8, width=6)

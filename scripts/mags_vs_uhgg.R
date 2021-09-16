# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mashdiff/")
dnadiff.uhgg = read.delim("mags_vs_uhgg-all_besthits.tsv", header=TRUE, stringsAsFactors = FALSE)
dnadiff.uhgg$ref = gsub(".fa", "", dnadiff.uhgg$ref)
checkm.data = read.csv("../qc_stats/all_checkm.csv")
checkm.data$genome = gsub(".fa", "", checkm.data$genome)
genome2species = read.delim("../../Project_resource/genomes-all_metadata.tsv", stringsAsFactors = FALSE)[,c("Genome","Species_rep")]
gtdb.tk = read.delim("../../Project_resource/species_taxonomy.tab")
colnames(gtdb.tk)[1] = "UHGG_species"

# prepare dataset
dnadiff.uhgg$UHGG_species = genome2species$Species_rep[match(dnadiff.uhgg$ref, genome2species$Genome)]
dnadiff.uhgg = merge(dnadiff.uhgg, gtdb.tk, by="UHGG_species")
dnadiff.uhgg[,c("UHGG_comp", "UHGG_cont")] = checkm.data[match(dnadiff.uhgg$UHGG_species, checkm.data$genome), 
                                                        c("completeness", "contamination")]
dnadiff.uhgg$UHGG_qs = dnadiff.uhgg$UHGG_comp - 5*dnadiff.uhgg$UHGG_cont
dnadiff.uhgg[,c("MAGs_comp", "MAGs_cont")] = checkm.data[match(dnadiff.uhgg$qry, checkm.data$genome), 
                                                         c("completeness", "contamination")]
dnadiff.uhgg$MAGs_qs = dnadiff.uhgg$MAGs_comp - 5*dnadiff.uhgg$MAGs_cont
dnadiff.uhgg$af_max = apply(dnadiff.uhgg[,c("ref_cov", "qry_cov")], 1, max)
dnadiff.uhgg$Match = ifelse(dnadiff.uhgg$af_max >= 80 & dnadiff.uhgg$ani >= 99, "Strain match",
                            ifelse(dnadiff.uhgg$af_max >= 30 & dnadiff.uhgg$ani >= 95, "Species match",
                                   "No match"))

# dnadiff scatterplot
dnadiff.plot = ggplot(dnadiff.uhgg, aes(x=qry_cov, y=ref_cov, colour=Match)) +
  geom_point(size=0.5, alpha=0.4) +
  theme_classic() +
  scale_colour_manual(values=c("tomato", "steelblue", "olivedrab3"), name="Match level") +
  ylab("% UHGG aligned") +
  xlab("% MAG aligned") +
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
  scale_fill_manual(values=rev(c("tomato", "steelblue", "olivedrab3")), name="Match level") +
  scale_y_continuous(breaks=c(0,4000,8000,12000,16000)) +
  ylab("Number of MAGs") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# arrange plot and save
ggarrange(ggarrange(dens.plot, dnadiff.plot, ncol=1, nrow=2, heights=c(1,4), align="v", legend = FALSE), 
          bar.plot, ncol=2, nrow=1, widths=c(2.5,1), common.legend = TRUE)
ggsave("dnadiff_vs_uhgg.pdf", height=6, width=8)

# save data
archaea.mags = dnadiff.uhgg[which(grepl("d__Archaea", dnadiff.uhgg$GTDB_lineage) & dnadiff.uhgg$Match != "No match"),"qry"]
unknown.mags = dnadiff.uhgg[which(dnadiff.uhgg$Match == "No match"),"qry"]
write.table(unknown.mags, file="unknown_mags.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(archaea.mags, file="archaea_mags.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
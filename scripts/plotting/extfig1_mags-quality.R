# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
checkm.stats = read.delim("prokaryotes/genomes_uhgg-latinbiota.tsv", stringsAsFactors = FALSE)
checkm.stats = checkm.stats[!grepl("GUT_GENOME", checkm.stats$Genome),]
checkm.stats$MIMAG = ifelse(checkm.stats$Completeness > 90 & checkm.stats$Contamination < 5, "Near complete", "Medium quality")
checkm.stats$MIMAG[which(checkm.stats$Completeness > 90 & checkm.stats$Contamination < 5 &
                           checkm.stats$rRNA_5S > 80 & checkm.stats$rRNA_16S > 80 & checkm.stats$rRNA_23S > 80 &
                           checkm.stats$tRNAs >= 18)] = "High quality"
checkm.stats$MIMAG = factor(checkm.stats$MIMAG, levels=c("Medium quality", "Near complete", "High quality"))

# plots
scatter.plot = ggplot(checkm.stats, aes(x=Completeness, y=Contamination, colour=MIMAG)) + 
  geom_point(size=0.8, alpha=0.5) + 
  scale_color_manual(values=c("steelblue", "darkolivegreen3", "darkgreen"), name="") + 
  guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8))) + 
  theme_bw() + 
  ylab("Contamination (%)\n") + 
  xlab("Completeness (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))

comp.dens = ggplot(checkm.stats, aes(x=Completeness, fill=MIMAG)) + 
  geom_density(alpha=0.8, colour=NA) + 
  geom_vline(xintercept = median(checkm.stats$Completeness), linetype="dashed") +
  scale_fill_manual(values=c("steelblue", "darkolivegreen3", "darkgreen")) + 
  theme_bw() + 
  guides(fill=FALSE) +
  ylab("\n") + 
  theme(panel.border = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_blank())

cont.dens = ggplot(checkm.stats, aes(x=Contamination, fill=MIMAG)) + 
  geom_density(alpha=0.8, colour=NA) + 
  geom_vline(xintercept = median(checkm.stats$Contamination), linetype="dashed") +
  scale_fill_manual(values=c("steelblue", "darkolivegreen3", "darkgreen")) + 
  theme_bw() + 
  coord_flip() +
  ylab("\n") +
  theme(panel.border = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.title.x = element_text(size=12)) + 
  theme(axis.text.y = element_blank())

# arrange plots
ggarrange(ggarrange(comp.dens, ncol=2, widths=c(4,1)), ggarrange(scatter.plot, cont.dens, ncol=2, nrow=1, widths=c(4,1), common.legend = TRUE, legend="bottom"), 
          ncol=1, nrow=2, heights=c(1,3.5))
ggsave("../../../Publications/2022-Latinbiota/Figures/extfig1/mag_quality.pdf", height=6, width=8)
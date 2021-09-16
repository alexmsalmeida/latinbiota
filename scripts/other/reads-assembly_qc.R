# load libraries
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/qc_stats/")
read.counts = read.delim("read_counts.tsv", stringsAsFactors = FALSE, header = FALSE)
assemb.stats = read.delim("assembly_stats.tsv", stringsAsFactors = FALSE, header = FALSE)

# plot stats
read.plot = ggplot(read.counts, aes(x=V2*2)) +
  geom_histogram(fill="lightgrey", colour="darkgrey") +
  theme_classic() +
  ylab("Number of datasets") +
  xlab("Number of reads") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

assemb.scatter = ggplot(assemb.stats, aes(x=V2, y=V3)) +
  geom_point(colour="darkgrey") +
  theme_classic() +
  ylab("Number of contigs") +
  xlab("Assembly length (bp)") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

assemb.dens = ggplot(assemb.stats, aes(x=V2)) +
  geom_density(fill="darkgrey", colour="white") +
  theme_classic() +
  ylab("Density") +
  xlab("Assembly length (bp)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank())

# arrange and save
ggarrange(read.plot, ggarrange(assemb.dens, assemb.scatter, ncol=1, nrow=2, heights=c(1,3), align = "v"), ncol=2, nrow=1)
ggsave(file="qc_stats.tiff", height=5, width=9, dpi=300) 
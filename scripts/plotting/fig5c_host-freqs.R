# load libraries
library(RColorBrewer)
library(ggplot2)
library(tidyr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
host.bacteria.tax.filt = read.delim("viruses/vcontact/virus_host_summary.tsv", header=TRUE, stringsAsFactors = FALSE)

# parse taxonomy
keep.taxa = sort(names(sort(table(host.bacteria.tax.filt$Family), decreasing=TRUE)[1:10]))
host.bacteria.tax.filt$Family = ifelse(host.bacteria.tax.filt$Family %in% keep.taxa, host.bacteria.tax.filt$Family, "Other")
order.taxa = c(names(sort(table(host.bacteria.tax.filt$Family[which(host.bacteria.tax.filt$Family != "Other")]), decreasing=TRUE)), "Other")
host.bacteria.tax.filt$Family = factor(host.bacteria.tax.filt$Family, levels=order.taxa)
taxa.colors = c(brewer.pal(10, "Set3"), "darkgrey")
names(taxa.colors) = order.taxa

# plot family counts
df = data.frame(table(host.bacteria.tax.filt$Family))
colnames(df) = c("Family", "Counts")
host.plot = ggplot(df, aes(x=Family, y=Counts, fill=Family)) +
  geom_bar(stat="identity", colour="darkgrey", alpha=0.8, size=0.3) +
  theme_classic() +
  coord_flip() +
  guides(fill="none") +
  ylab("Number of viral species") +
  scale_fill_manual(values=taxa.colors) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank())
ggsave("../../../Publications/2022-Latinbiota/Figures/figure5/host_counts.pdf", width=7, height=5)

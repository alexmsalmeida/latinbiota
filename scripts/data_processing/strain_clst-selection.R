# load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

# set wkdir
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")

# load data
mash.file = read.csv("prokaryotes/oscillospiraceae/Cdb.csv", stringsAsFactors = FALSE)
mash.file$Source = NA
mash.file$Source[grep("bin", mash.file$genome)] = "Latinbiota"
mash.file$Source[which(is.na(mash.file$Source))] = "UHGG"
mash.file.latin = mash.file[grep("bin", mash.file$genome),]

# group by cluster and metadata
mash.grouped = mash.file %>%
  group_by(Source, secondary_cluster) %>%
  tally() %>%
  as.data.frame()

# calculate proportions
for (n in 1:nrow(mash.grouped)) {
  cluster = mash.grouped$secondary_cluster[n]
  mash.grouped$prop[n] = mash.grouped[n,"n"]/sum(mash.grouped[which(mash.grouped$secondary_cluster == cluster),"n"])*100
}
mash.grouped = mash.grouped[order(mash.grouped$prop, mash.grouped$n, decreasing=TRUE),]

# plot proportions of all
top.clst = names(sort(table(mash.file.latin$secondary_cluster), decreasing=TRUE)[1:20])
bar.plot = ggplot(mash.grouped, aes(x=secondary_cluster, y=n, fill=Source)) +
  geom_bar(stat="identity", alpha=0.8) +
  scale_x_discrete(limits=top.clst) +
  scale_fill_manual(values=c("darkred", "steelblue")) + 
  ylab("Number of genomes") + 
  xlab("Cluster (90% identity)") +
  ggtitle("Clusters ordered by total number of Latinbiota genomes") +
  theme_classic() + 
  theme(legend.position="right", legend.text = element_text(size=12), legend.title = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# plot proportions of latinbiota
latin.clst = mash.grouped[which(mash.grouped$Source == "Latinbiota" & mash.grouped$n > 10)[1:20],"secondary_cluster"]
latin.plot = ggplot(mash.grouped, aes(x=secondary_cluster, y=prop, fill=Source)) +
  geom_bar(stat="identity", alpha=0.8) +
  geom_text(aes(label=n), position = position_stack(vjust=0.5)) +
  scale_x_discrete(limits=latin.clst) +
  scale_fill_manual(values=c("darkred", "steelblue")) + 
  ylab("Proportion of genomes (%)") + 
  xlab("Cluster (90% identity)") +
  ggtitle("Clusters ordered by proportion of Latinbiota genomes") +
  theme_classic() + 
  theme(legend.position="right", legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# arrange and save
ggarrange(bar.plot, latin.plot, nrow=2, common.legend = TRUE, legend = "right")

# save selected clusters to file
selected.clst = names(sort(table(mash.file.latin$secondary_cluster), decreasing=TRUE)[1:5])
write.table(selected.clst, file="prokaryotes/oscillospiraceae/selected_clsts.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

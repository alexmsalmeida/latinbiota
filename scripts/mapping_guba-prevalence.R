# load libraries
library(tidyr)
library(reshape2)
library(data.table)
library(phytools)

# function
group_by_metadata = function(x, y) {
  if(y == "Country"){
    bwa.group = aggregate(. ~ Country, data=x[,-which(colnames(x) == "Continent")], FUN=sum)
  } else {
    bwa.group = aggregate(. ~ Continent, data=x[,-which(colnames(x) == "Country")], FUN=sum)
  }
  meta.counts = data.frame(table(metadata[which(metadata$Run %in% rownames(x)),y]))
  location.prop = data.frame(t(bwa.group[,-1]/meta.counts$Freq*100))
  colnames(location.prop) = meta.counts$Var1
  return(location.prop)
}

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/phylogeny/gubaphage/")
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
latinbiota = sort(unique(metadata[which(metadata$Study == "Latinbiota"),"Country"]))
latinbiota.samples = metadata[which(metadata$Study == "Latinbiota"),"Run"]
guba.genomes = read.tree("gubaphage_tree.nwk")$tip.label
bwa.data = as.data.frame(fread("../../mapping/gpd-combined/summary_tables/viruses-above0_binary_world-adult-healthy.csv"))
rownames(bwa.data) = bwa.data$Sample
bwa.data = bwa.data[,-c(1,ncol(bwa.data))]
bwa.data = bwa.data[,which(colnames(bwa.data) %in% guba.genomes)]

# calculate total freqs and stats
total.freqs = data.frame(Freq=colSums(bwa.data))
global.prop = length(which(rowSums(bwa.data) > 0))/nrow(bwa.data)*100
latinbiota.prop = length(which(rowSums(bwa.data[which(rownames(bwa.data) %in% latinbiota.samples),]) > 0))/nrow(bwa.data[which(rownames(bwa.data) %in% latinbiota.samples),])*100
n_world.guba = length(which(rowSums(bwa.data) > 0))
n_world.noguba = nrow(bwa.data) - n_world.guba
n_latinbiota.guba = length(which(rowSums(bwa.data[which(rownames(bwa.data) %in% latinbiota.samples),]) > 0))
n_latinbiota.noguba = nrow(bwa.data[which(rownames(bwa.data) %in% latinbiota.samples),]) - n_latinbiota.guba
chisq.test(rbind(c(n_world.guba, n_world.noguba), c(n_latinbiota.guba, n_latinbiota.noguba)))

# group by metadata
bwa.data$Country = metadata[rownames(bwa.data), "Country"]
bwa.data$Continent = metadata[rownames(bwa.data), "Continent"]
bwa.continent = group_by_metadata(bwa.data, "Continent")
bwa.countries = group_by_metadata(bwa.data, "Country")
bwa.latinbiota = bwa.countries[,latinbiota]

# save data
write.table(cbind(rownames(total.freqs), total.freqs), file="guba_freqs.txt", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")
write.table(cbind(rownames(bwa.latinbiota), bwa.latinbiota), file="latinbiota_props.txt", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")
write.table(cbind(rownames(bwa.continent), bwa.continent), file="continent_props.txt", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")

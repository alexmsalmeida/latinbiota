# load libraries
library(data.table)
library(matrixStats)

# function
read_bwa = function (x) {
  bwa.data = as.data.frame(fread(x, check.names = FALSE))
  rownames(bwa.data) = bwa.data$Genome
  bwa.data = bwa.data[,metadata$Run]
  return(bwa.data)
}

# set directory
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")

# define metadata criteria
metadata = read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
metadata = metadata[which(metadata$`Health state` == "Healthy" & metadata$`Age group` == "Adult" & metadata$Antibiotics != "Yes"),]; variable = "Continent"
#metadata = metadata[which(metadata$Study == "Latinbiota"),]; variable = "Country"

# load mapping data
bwa.cov.est = read_bwa("prokaryotes/mapping/raw/bwa_cov-est.csv")
bwa.cov.exp = read_bwa("prokaryotes/mapping/raw/bwa_cov-exp.csv")

# define prevalence cut-offs
cov.est = 10 # Percentage
cov.ratio = 0.3 # Ratio
prop_samples = 0 # Ratio

bwa.prev = bwa.cov.est
bwa.prev[bwa.cov.est < cov.est | bwa.cov.est/bwa.cov.exp < cov.ratio] = 0
bwa.prev[bwa.prev != 0] = 1
species = names(which(rowSums(bwa.prev) > ncol(bwa.prev)*prop_samples))

# subset by species and samples
bwa.binary = bwa.prev[species,]
bwa.binary = data.frame(t(bwa.binary[,names(which(colSums(bwa.binary) > 0))]), check.names=FALSE)

# transform count data
bwa.counts = data.frame(t(read_bwa("prokaryotes/mapping/raw/bwa_counts-unique.csv")), check.names = FALSE)
bwa.counts = bwa.counts[rownames(bwa.binary), colnames(bwa.binary)]
bwa.counts[bwa.binary == 0] = 0
bwa.binary = bwa.counts

# add metadata variable
meta.var = make.names(metadata[rownames(bwa.binary), variable])
bwa.binary$Variable = meta.var

# save final table
bwa.binary = data.frame(cbind(rownames(bwa.binary), bwa.binary), check.names = FALSE)
colnames(bwa.binary)[1] = "Sample"
write.table(bwa.binary, file="prokaryotes/mapping/summary_tables/euproks-above0_counts-transf_all.csv", sep=",", quote=FALSE, row.names=FALSE)

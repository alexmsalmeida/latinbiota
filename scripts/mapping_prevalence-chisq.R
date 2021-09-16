# load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(effectsize)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/uhgg-combined")
bwa.data = as.data.frame(read.csv("summary_tables/euproks-above0_binary_world-adult-healthy.csv", 
                                  row.names=1, check.names=FALSE, stringsAsFactors = FALSE))
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
rownames(metadata) = metadata$Run
bwa.data$Variable = ifelse(metadata[rownames(bwa.data),"Study"] == "Latinbiota", "Latinbiota", "Public")

# aggregate data
bwa.agg = aggregate(. ~ Variable, data=bwa.data, FUN=sum)
bwa.prev = as.data.frame(t(bwa.agg[,-1]))
bwa.miss = as.data.frame(t(abs(t(bwa.prev)-as.vector(table(bwa.data$Variable)))))
bwa.merge = cbind(bwa.prev, bwa.miss)
colnames(bwa.merge) = c("Latinbiota_yes", "Public_yes", "Latinbiota_no", "Public_no")

# perform chisq test
bwa.merge$effsize = apply(bwa.merge, 1, function(x) cramers_v(matrix(x,nrow = 2))$Cramers_v)
bwa.merge$p_value = apply(bwa.merge[,-ncol(bwa.merge)], 1, function(x) chisq.test(matrix(x,nrow = 2))$p.value)
bwa.merge$fdr = p.adjust(bwa.merge$p_value, method="fdr")
bwa.sign = bwa.merge[which(bwa.merge$fdr < 0.05),]
bwa.sign$result = ifelse(bwa.sign[,1]/sum(bwa.sign[,c(1,3)]) > bwa.sign[,2]/sum(bwa.sign[,c(2,4)]), "Higher", "Lower")

# confirm significance by permutation
n_permutations = 1000
permute_res = function(x){
  cat("Iteration", x, "\n")
  bwa.subs = bwa.data[,c("Variable", rownames(bwa.sign))]
  bwa.perm = rbind(bwa.subs[which(bwa.subs$Variable == "Latinbiota"),], 
                   bwa.subs[sample(which(bwa.subs$Variable == "Public"), length(which(bwa.subs$Variable == "Latinbiota"))),])
  bwa.perm.agg = aggregate(. ~ Variable, data=bwa.perm, FUN=sum)
  bwa.perm.prev = as.data.frame(t(bwa.perm.agg[,-1]))
  bwa.perm.miss = as.data.frame(t(abs(t(bwa.perm.prev)-as.vector(table(bwa.perm$Variable)))))
  bwa.perm.merge = cbind(bwa.perm.prev, bwa.perm.miss)
  colnames(bwa.perm.merge) = c("Latinbiota_yes", "Public_yes", "Latinbiota_no", "Public_no")
  perm.result = apply(bwa.perm.merge,1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
  perm.result = p.adjust(perm.result, method="fdr")
  names(perm.result) = rownames(bwa.perm.merge)
  return(perm.result)
}

set.seed(1) 
result = sapply(1:n_permutations, permute_res)

# combine final results
bwa.sign$permutation = rowSums(result >= 0.05)/n_permutations
bwa.fi = bwa.sign[which(bwa.sign$permutation < 0.05),]
bwa.fi = cbind(rownames(bwa.fi), bwa.fi)
colnames(bwa.fi)[1] = "Genome"
write.table(bwa.fi, file="chisq_euproks_latinbiota.tsv", sep="\t", quote=FALSE, row.names=FALSE)

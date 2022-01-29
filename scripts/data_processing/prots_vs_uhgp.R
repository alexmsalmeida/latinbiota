# load libraries
library(ggplot2)
library(ggrastr)
library(ggpubr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
diamond.res = read.delim("prokaryotes/functions/latinbiota_vs_uhgp.tsv", header=FALSE, stringsAsFactors = FALSE)
colnames(diamond.res) = c("qry", "ref", "pid", "qlen", "qstart", "qend", "slen", "sstart", "send", "bitscore", "evalue")
diamond.res$qry_cov = (diamond.res$qend-diamond.res$qstart)/diamond.res$qlen*100
diamond.res$ref_cov = (diamond.res$send-diamond.res$sstart)/diamond.res$slen*100
diamond.res$cov_max = apply(diamond.res[,c("ref_cov", "qry_cov")], 1, max)
diamond.res$match = ifelse(diamond.res$cov_max >= 80 & diamond.res$pid >= 90, "Known protein", "Novel protein")
write.table(diamond.res$qry[which(diamond.res$match == "Novel protein")], file="prokaryotes/functions/novel_proteins.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)

# load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")
novel.proteins = scan("prokaryotes/functions/novel_proteins.txt", what="")
file = "prokaryotes/functions/annotations/CAZy_annotations.tsv"
no_col = max(count.fields(file, sep = "\t"))
cazy.data = read.table(file, sep="\t", fill=TRUE, header = F, col.names=c(1:no_col))
cazy.descr = read.delim("prokaryotes/functions/cazy_descriptions.tsv", header = FALSE, stringsAsFactors = FALSE)
cazy.descr$V2 = ifelse(str_count(cazy.descr$V2, pattern="/") > 1, "Unspecific", cazy.descr$V2)

# extract KOs
cazy.novel = as.vector(as.matrix(cazy.data[which(cazy.data$X1 %in% novel.proteins),-1]))
cazy.novel = cazy.novel[which(cazy.novel != "")]
cazy.novel.df = data.frame(sort(table(cazy.novel), decreasing=TRUE))
colnames(cazy.novel.df) = c("CAZy", "Novel_present")

cazy.known = as.vector(as.matrix(cazy.data[which(!cazy.data$X1 %in% novel.proteins),-1]))
cazy.known = cazy.known[which(cazy.known != "")]
cazy.known.df = data.frame(sort(table(cazy.known), decreasing=TRUE))
colnames(cazy.known.df) = c("CAZy", "Known_present")

cazy.all = merge(cazy.novel.df, cazy.known.df, by="CAZy", all = TRUE)
cazy.all[is.na(cazy.all)] = 0

cazy.all$Novel_absent = sum(cazy.all$Novel_present)-cazy.all$Novel_present
cazy.all$Known_absent = sum(cazy.all$Known_present)-cazy.all$Known_present

# calculate stats
cols = c("Novel_present", "Known_present", "Novel_absent", "Known_absent")
cazy.all$pvalue = apply(cazy.all[,cols],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
cazy.all$fdr = p.adjust(cazy.all$pvalue, method="fdr")
cazy.all$result = ifelse(cazy.all[,2]/sum(cazy.all[,c(2,4)]) > cazy.all[,3]/sum(cazy.all[,c(3,5)]), "Higher", "Lower")
cazy.sign = cazy.all[which(cazy.all$fdr < 0.05),]
cazy.sign = cazy.sign[order(cazy.sign$fdr, decreasing=FALSE),]
cazy.sign$funcs = cazy.descr[match(cazy.sign$CAZy, cazy.descr$V1),"V2"]
cazy.sign$funcs = ifelse(is.na(cazy.sign$funcs), paste(cazy.sign$CAZy, "(Unspecific)"), 
                         paste(cazy.sign$CAZy, " (", cazy.sign$funcs, ")", sep=""))

# plot significant functions
cazy.plot = ggplot(cazy.sign, aes(x=CAZy, y=Novel_present, fill=result)) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_flip() +
  theme_classic() +
  ylab("Number of new proteins") +
  scale_x_discrete(limits=cazy.sign[order(cazy.sign$Novel_present, decreasing=TRUE)[1:23],"CAZy"],
                   labels=cazy.sign[order(cazy.sign$Novel_present, decreasing=TRUE)[1:23],"funcs"]) +
  scale_fill_manual(values=c("steelblue", "tomato"), labels=c("Enriched in novel", "Enriched in known"), name="") +
  theme(legend.position="bottom") + 
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12))
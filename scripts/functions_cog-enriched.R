# load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/functions/")
novel.proteins = scan("novel_proteins.txt", what="")
cogs = read.delim("eggNOG_annotations.tsv", header=FALSE, stringsAsFactors = FALSE)[,c(1,7)]
cog.descr = read.delim("cog_descriptions.tsv", stringsAsFactors = FALSE)

# extract COGs
cog.novel = unlist(strsplit(cogs$V7[which(cogs$V1 %in% novel.proteins)], split=""))
cog.novel = cog.novel[which(cog.novel != "-")]
cog.novel.df = data.frame(sort(table(cog.novel), decreasing=TRUE))
colnames(cog.novel.df) = c("COG", "Novel_present")

cog.known = unlist(strsplit(cogs$V7[which(!cogs$V1 %in% novel.proteins)], split=""))
cog.known = cog.known[which(cog.known != "-")]
cog.known.df = data.frame(sort(table(cog.known), decreasing=TRUE))
colnames(cog.known.df) = c("COG", "Known_present")

cog.all = merge(cog.novel.df, cog.known.df, by="COG", all = TRUE)
cog.all[is.na(cog.all)] = 0

cog.all$Novel_absent = sum(cog.all$Novel_present)-cog.all$Novel_present
cog.all$Known_absent = sum(cog.all$Known_present)-cog.all$Known_present

# calculate stats
cols = c("Novel_present", "Known_present", "Novel_absent", "Known_absent")
cog.all$pvalue = apply(cog.all[,cols],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
cog.all$fdr = p.adjust(cog.all$pvalue, method="fdr")
cog.all$result = ifelse(cog.all[,2]/sum(cog.all[,c(2,4)]) > cog.all[,3]/sum(cog.all[,c(3,5)]), "Higher", "Lower")
cog.sign = cog.all[which(cog.all$fdr < 0.05),]
cog.sign = cog.sign[order(cog.sign$fdr, decreasing=FALSE),]
cog.sign$subcat = cog.descr[match(cog.sign$COG, cog.descr$subcat_id),"subcat_name"]
cog.sign$cat = cog.descr[match(cog.sign$COG, cog.descr$subcat_id),"cat_name"]

# plot significant functions
cog.plot = ggplot(cog.sign, aes(x=subcat, y=Novel_present, fill=result)) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_flip() +
  theme_classic() +
  ylab("Number of new proteins") +
  scale_x_discrete(limits=cog.sign[order(cog.sign$Novel_present, decreasing=TRUE),"subcat"]) +
  scale_fill_manual(values=c("steelblue", "tomato"), labels=c("Enriched in novel", "Enriched in known"), name="") +
  theme(legend.position="bottom") + 
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12))
ggsave("cog_enriched.png", dpi=300, height=5, width=8)
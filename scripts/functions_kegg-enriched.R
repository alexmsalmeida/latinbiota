# load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/functions/")
novel.proteins = scan("novel_proteins.txt", what="")
kegg.ontology = read.delim("kegg_ontology/kegg_classes.tsv", header=FALSE, stringsAsFactors = FALSE, sep="\t")
kegg.ontology = separate(kegg.ontology, V4, into = c("KO", "annotation"), sep = "\\s", extra = "merge")
kegg.ontology = separate(kegg.ontology, annotation, into = c("Gene", "Annotation"), sep = "; ", extra = "merge")
kegg.ontology$Annotation = gsub(" \\[.*", "", kegg.ontology$Annotation)
file = "KOFam_annotations.tsv"
no_col = max(count.fields(file, sep = "\t"))
kegg.data = read.table(file, sep="\t", fill=TRUE, header = F, col.names=c(1:no_col))

# extract KOs
kegg.novel = as.vector(as.matrix(kegg.data[which(kegg.data$X1 %in% novel.proteins),-1]))
kegg.novel = kegg.novel[which(kegg.novel != "")]
kegg.novel.df = data.frame(sort(table(kegg.novel), decreasing=TRUE))
colnames(kegg.novel.df) = c("KO", "Novel_present")

kegg.known = as.vector(as.matrix(kegg.data[which(!kegg.data$X1 %in% novel.proteins),-1]))
kegg.known = kegg.known[which(kegg.known != "")]
kegg.known.df = data.frame(sort(table(kegg.known), decreasing=TRUE))
colnames(kegg.known.df) = c("KO", "Known_present")

kegg.all = merge(kegg.novel.df, kegg.known.df, by="KO", all = TRUE)
kegg.all[is.na(kegg.all)] = 0

kegg.all$Novel_absent = sum(kegg.all$Novel_present)-kegg.all$Novel_present
kegg.all$Known_absent = sum(kegg.all$Known_present)-kegg.all$Known_present

# calculate stats
cols = c("Novel_present", "Known_present", "Novel_absent", "Known_absent")
kegg.all$pvalue = apply(kegg.all[,cols],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
kegg.all$fdr = p.adjust(kegg.all$pvalue, method="fdr")
kegg.all$result = ifelse(kegg.all[,2]/sum(kegg.all[,c(2,4)]) > kegg.all[,3]/sum(kegg.all[,c(3,5)]), "Higher", "Lower")
kegg.sign = kegg.all[which(kegg.all$fdr < 0.05),]
kegg.sign = kegg.sign[order(kegg.sign$fdr, decreasing=FALSE),]
kegg.sign$Annotation = kegg.ontology[match(kegg.sign$KO, kegg.ontology$KO),"Annotation"]

# plot significant functions
kegg.plot = ggplot(kegg.sign, aes(x=KO, y=Novel_present, fill=result)) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_flip() +
  theme_classic() +
  ylab("Number of new proteins") +
  scale_x_discrete(limits=kegg.sign[order(kegg.sign$Novel_present, decreasing=TRUE)[1:25],"KO"],
                   labels=kegg.sign[order(kegg.sign$Novel_present, decreasing=TRUE)[1:25],"Annotation"]) +
  scale_fill_manual(values=c("steelblue", "tomato"), labels=c("Enriched in novel", "Enriched in known"), name="") +
  theme(legend.position="bottom") + 
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=12))
ggsave("kegg_enriched.png", dpi=300, height=5, width=6)
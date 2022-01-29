# load libraries
library(data.table)
library(ggpubr)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota")
gwas.results = read.csv("prokaryotes/oscillospiraceae/sign-genes.tsv", stringsAsFactors = FALSE)
gwas.latinbiota = gwas.results[which(gwas.results[,4]/gwas.results[,6] > gwas.results[,5]/gwas.results[,7] & gwas.results[,13] < 0.05),]
gwas.uhgg = gwas.results[which(gwas.results[,4]/gwas.results[,6] < gwas.results[,5]/gwas.results[,7] & gwas.results[,13] < 0.05),]
sign.proteins = read.table("prokaryotes/oscillospiraceae/annotations_sign-genes.tsv", header=FALSE, stringsAsFactors = FALSE, sep = "\t")[,c(1,7)]
sign.proteins = sign.proteins[which(sign.proteins$V1 %in% gwas.latinbiota$Gene),]
pan.genome = read.table("prokaryotes/oscillospiraceae/annotations_pan-genome.tsv", header=FALSE, stringsAsFactors = FALSE, sep = "\t")[,c(1,7)]
cog.cats = read.delim("prokaryotes/functions/cog_descriptions.tsv", stringsAsFactors = FALSE)

# extract COGs
cog.sign = unlist(strsplit(sign.proteins$V7, split=""))
cog.sign = cog.sign[which(cog.sign != "-")]
cog.sign.df = data.frame(sort(table(cog.sign), decreasing=TRUE))
colnames(cog.sign.df) = c("subcat_id", "Sign_present")

cog.all = unlist(strsplit(pan.genome$V7, split=""))
cog.all = cog.all[which(cog.all != "-")]
cog.all.df = data.frame(sort(table(cog.all), decreasing=TRUE))
colnames(cog.all.df) = c("subcat_id", "Total_present")

cog.all = merge(cog.sign.df, cog.all.df, by="subcat_id", all = TRUE)
cog.all[is.na(cog.all)] = 0

cog.all$Sign_absent = sum(cog.all$Sign_present)-cog.all$Sign_present
cog.all$Total_absent = sum(cog.all$Total_present)-cog.all$Total_present

# calculate stats
cols = c("Sign_present", "Total_present", "Sign_absent", "Total_absent")
cog.all$pvalue = apply(cog.all[,cols],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
cog.all$fdr = p.adjust(cog.all$pvalue, method="fdr")
cog.all$result = ifelse(cog.all[,2]/sum(cog.all[,c(2,4)]) > cog.all[,3]/sum(cog.all[,c(3,5)]), "Higher", "Lower")

# add COG description
cog.df = merge(cog.all, cog.cats, by="subcat_id")
cog.df = cog.df[which(cog.df$cat_name != "Poorly characterized" & cog.df$Sign_present > 0),]

# plot
bar.plot = ggplot(cog.df, aes(y=Sign_present, x=reorder(subcat_name, -Sign_present), fill=cat_name)) + 
  geom_bar(stat="identity", alpha=0.8, size=0.1, colour="darkgrey", width=0.7) + 
  theme_classic() + 
  coord_flip() + 
  scale_fill_manual(values=c("tomato", "steelblue", "darkgreen"), name="Functional category") +
  ylab("Number of genes") +
  theme(legend.position = "top", legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=14)) +
  theme(plot.margin=unit(c(1,1,0.5,1.2),"cm")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank())
bar.nolegend = bar.plot + theme(legend.position = "none")
legend = get_legend(bar.plot)
ggarrange(bar.nolegend, legend, nrow=2, heights=c(17,1))
ggsave(filename="../../../Publications/2022-Latinbiota/Figures/extfig9/cog_enriched-latinbiota.pdf", width=12, height=7)

# load libraries
library(pheatmap)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/machine_learning/")
model = readRDS("latinbiota_binary/glmnet_109_model.Rds")
taxa.all = read.delim("../../viral_prediction/DemoVir_assignments.txt", row.names=1)
metadata = read.delim("../../metadata/metadata_complete.tsv", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
bwa.data = as.data.frame(read.csv("../viruses-above10_binary_latinbiota-all.csv", row.names=1, check.names=FALSE, stringsAsFactors = FALSE))

# subset features
#best.features = as.vector(model$feature_importance$names[which(model$feature_importance$perf_metric_diff > 0)])
best.features = as.vector(top_n(model$feature_importance, n=100, wt=perf_metric_diff)$names)
best.features = gsub("_1$", "", best.features)
feat.heat = bwa.data[,best.features]
feat.heat = bwa.data[,names(sort(colSums(feat.heat), decreasing=TRUE))]

# prepare annotation
taxa.all = taxa.all[which(rownames(taxa.all) %in% colnames(feat.heat)),]
feat.taxa = data.frame(Family=taxa.all[,"Family"], row.names=rownames(taxa.all))
taxa.colors = brewer.pal(length(unique(feat.taxa$Family)), "Set3")
names(taxa.colors) = unique(feat.taxa$Family)
metadata = metadata[which(metadata$Run %in% rownames(feat.heat)),]
samp.geo = data.frame(Location=metadata[,"Country"], row.names=metadata[,"Run"])
geo.colors =c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#b2450a", "#278d29", "#3b35ce")
names(geo.colors) = sort(as.vector(unique(samp.geo$Location)))
annot.colors = list(Family = taxa.colors, Location=geo.colors)

# plot heatmap
pheatmap(feat.heat, annotation_col = feat.taxa, annotation_row = samp.geo, annotation_colors = annot.colors,
         cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = c("grey90", "steelblue"))

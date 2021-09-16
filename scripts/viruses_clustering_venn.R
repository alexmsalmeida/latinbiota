# load library
library(VennDiagram)
library(ggplot2)

# gpd and latinbiota
latinbiota = 15881
gpd = 39032
gpd.combined = 48130
overlap = gpd + latinbiota - gpd.combined
latinbiota_increase = (gpd.combined - gpd)/gpd*100

# double venn
pdf(file="~/Documents/ESPOD/Analyses/Project_latinbiota/viral_prediction/gpd-combined_venn.pdf", height=5, width=6)
draw.pairwise.venn(area1 = gpd, area2 = latinbiota, cross.area=overlap,
                   category = c("GPD", "Latinbiota"),
                   lty = rep("blank", 2), fill = c("plum3", "darkorange"), 
                   alpha = rep(0.4, 2),
                   cex=1.5, fontfamily = rep("sans", 3),
                   cat.cex = 1.7, cat.fontfamily = rep("sans", 2), 
                   euler.d = TRUE, scaled = TRUE)
dev.off()
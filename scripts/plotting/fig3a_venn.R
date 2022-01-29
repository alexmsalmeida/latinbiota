# load library
library(VennDiagram)
library(ggplot2)

# setwd
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/")

# gpd and latinbiota
latinbiota = 17167
gpd = 39055
gpd.combined = 49026
overlap = gpd + latinbiota - gpd.combined
latinbiota_increase = (gpd.combined - gpd)/gpd*100

# double venn
pdf(file="../../../Publications/2022-Latinbiota/Figures/figure3/venn_diagram.pdf", height=5, width=6)
draw.pairwise.venn(area1 = gpd, area2 = latinbiota, cross.area=overlap,
                   category = c("GPD", "Latinbiota"),
                   lty = rep("blank", 2), fill = c("plum3", "darkorange"), 
                   alpha = rep(0.4, 2),
                   cex=1.5, fontfamily = rep("sans", 3),
                   cat.cex = 1.7, cat.fontfamily = rep("sans", 2), 
                   euler.d = TRUE, scaled = TRUE)
dev.off()
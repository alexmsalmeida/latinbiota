# load library
library(curatedMetagenomicData)
library(readxl)

# load data
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/metadata")

# my metadata
my.meta = as.data.frame(read_excel("metadata_complete.xlsx"))
rownames(my.meta) = my.meta$Run
my.meta.healthy = my.meta[which(my.meta$`Health state` == "Healthy" & my.meta$`Age group` == "Adult" & my.meta$Antibiotics != "Yes"),"Run"]

# curatedMetagenomicData
all.meta = combined_metadata
all.meta = all.meta[which(all.meta$NCBI_accession %in% my.meta$Run),]
all.meta.healthy = all.meta[which(all.meta$disease == "healthy" & all.meta$age_category == "adult" & all.meta$antibiotics_current_use != "yes"),"NCBI_accession"]

# compare datasets
all.meta.subs = all.meta[which(all.meta$NCBI_accession %in% my.meta.healthy),]
to_exclude = all.meta.subs[which(all.meta.subs$disease != "healthy" & !is.na(all.meta.subs$disease) 
                            | all.meta.subs$age_category != "adult" & !is.na(all.meta.subs$age_category) & all.meta.subs$age < 18 | all.meta.subs$antibiotics_current_use == "yes"),"NCBI_accession"]

my.meta.subs = my.meta[which(my.meta$Run %in% all.meta.healthy),]
to_add = my.meta.subs[which(my.meta.subs$`Age group` == "NA" | my.meta.subs$`Health state` == "NA"),"Run"]

# reclassify data
my.meta[to_exclude, "Health state"] = "Diseased"
my.meta[to_add, "Health state"] = "Healthy"
my.meta[to_add, "Age group"] = "Adult"
write.table(my.meta, file="metadata_complete_revised.tsv", sep="\t", quote=FALSE, row.names=FALSE)


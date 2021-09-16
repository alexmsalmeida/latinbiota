# Latinbiota - gut microbiome diversity in Latin America
Central repository with scripts and workflows for the analysis of the gut microbiome of the Latinbiota cohort

Data available in: http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/latinbiota/

## Main analysis workflows

[VirSearch](https://github.com/alexmsalmeida/virsearch) - Searching viral sequences in metagenomes

[GubaScreen](https://github.com/alexmsalmeida/gubascreen) - Protocol for detecting Gubaphage genomes among a set of predicted viral sequences

[MAGscreen](https://github.com/alexmsalmeida/magscreen) - Workflow for the identification and curation of new microbial species

[metaMAP](https://github.com/alexmsalmeida/metamap) - Mapping sequencing reads to reference genomes for detecting microbial species

[ml-microbiome](https://github.com/alexmsalmeida/ml-microbiome) - Machine learning pipeline for classification of microbiome data.

## Statistics and visualization

All of the following scripts are available in the `scripts/` folder of this repo.

### Sample distribution and metadata:

* `samples_map_all.R`
* `samples_map_latinbiota.R`
* `samples_metadata_latinbiota.R`

### MAGs quality control and comparison with [UHGG](https://www.nature.com/articles/s41587-020-0603-3):

* `mags_checkm-stats.R`
* `mags_vs_uhgg.R`
* `mags_qs-comparison.R`
* `bacteria_phylo-diversity.R`

### Viral prediction, clustering and classification:

* `viruses_predicted_upset.R`
* `viruses_taxonomy.R`
* `viruses_cov-threshold.R`
* `viruses_clustering_venn.R`

### Diversity analysis

* `mapping_cutoffs.R`
* `mapping_kingdom-classification.R`
* `mapping_alpha-div.R`
* `mapping_euks-prevalence.R`
* `mapping_feat-heatmap_gpd.R`
* `mapping_feat-heatmap_uhgg.R`
* `mapping_guba-prevalence.R`
* `mapping_ml-summary.R`
* `mapping_pca.R`
* `mapping_prevalence-chisq.R`
* `mapping_prevalence-upset_gpd.R`
* `mapping_prevalence-upset_uhgg.R`


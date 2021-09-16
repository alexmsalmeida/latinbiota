# load libraries
library(mikropml)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# function
get_sens_spec_lookup = function(country, data){
  selected = data %>%
    select(!!as.name(country), observed)
  
  total = selected %>%
    count(observed) %>%
    pivot_wider(names_from=observed, values_from=n)
  
  selected %>%
    arrange(desc(!!as.name(country))) %>%
    mutate(is_country = (observed == country),
           tp = cumsum(is_country),
           fp = cumsum(!is_country),
           sensitivity = tp / as.numeric(total[country]),
           fpr = fp / (sum(total)-as.numeric(total[country])),
           specificity = 1-fpr) %>%
    add_column(class = country) %>%
    select(sensitivity, specificity, fpr, class)
}

# load all models
setwd("~/Documents/ESPOD/Analyses/Project_latinbiota/mapping/")
model.euproks.world = read.csv("uhgg-combined/machine_learning/euproks-above0_binary_world-adult-healthy.csv", stringsAsFactors = FALSE)
model.euproks.world$Database = "UHGG"
model.euproks.world$Samples = "Global"
model.viruses.world = read.csv("gpd-combined/machine_learning/viruses-above0_binary_world-adult-healthy.csv", stringsAsFactors = FALSE)
model.viruses.world$Database = "GPD"
model.viruses.world$Samples = "Global"
model.euproks.latin = read.csv("uhgg-combined/machine_learning/euproks-above0_binary_latinbiota-all.csv", stringsAsFactors = FALSE)
model.euproks.latin$Database = "UHGG"
model.euproks.latin$Samples = "Latinbiota"
model.viruses.latin = read.csv("gpd-combined/machine_learning/viruses-above0_binary_latinbiota-all.csv", stringsAsFactors = FALSE)
model.viruses.latin$Database = "GPD"
model.viruses.latin$Samples = "Latinbiota"
model.data = rbind(model.euproks.world, model.viruses.world, model.euproks.latin, model.viruses.latin)
model.data = model.data[which(model.data$Samples == "Global"),]

# plot score distributions
ml.boxplots = ggplot(model.data, aes(x=method, y=Kappa, fill=Database)) +
  geom_point(alpha=0.6, position=position_jitterdodge(jitter.width=0.05, dodge.width=0.5), size=0.5, 
                    colour="darkgrey") +
  geom_boxplot(alpha=0.5, outlier.colour = NA, colour="black", width=0.5) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(limits=c("svmRadial", "rf", "glmnet"), labels=c("SVM", "Random\nForest", "Logistic\nRegression")) +
  scale_fill_manual(values=c("steelblue", "darkgreen")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  theme(legend.position="top", legend.box = "horizontal", legend.text=element_text(size=12),
        legend.title = element_text(size=12, face="bold")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))
ggsave(file="ml_boxplot_latinbiota.png", dpi=300, height=5, width=6)

# select best model
best.model = top_n(model.data, n=1, wt=Kappa)
model = readRDS("uhgg-combined/machine_learning/euproks-above0_binary_world-adult-healthy_glmnet-106.Rds"); variable = "Continent"
#model = readRDS("gpd-combined/machine_learning/viruses-above0_binary_latinbiota-all_glmnet-105.Rds"); variable = "Country"

# prepare dataframe for ROC
prob = predict(model$trained_model, model$test_data, type="prob")
observed = model$test_data$Variable
prob_obs = bind_cols(prob, observed=observed)
roc.data = map_dfr(.x=unique(prob_obs$observed), .f=get_sens_spec_lookup, prob_obs)

# plot ROC curve
roc.curve = ggplot(roc.data, aes(x=1-specificity, y=sensitivity, colour=class)) +
  geom_line() +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  ylim(0,1) +
  xlim(0,1) +
  scale_colour_manual(values=c("#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", 
                               "#FC8D62", "#b2450a", "#3b35ce", "darkgrey"),name=variable) +
  ylab("True Positive Rate") + 
  xlab("False Positive Rate") +
  theme(legend.position="top", legend.box = "horizontal", legend.text=element_text(size=12),
        legend.title = element_text(size=12, face="bold")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))

# arrange and plot
ggarrange(ml.boxplots, roc.curve)
ggsave(file="ml_results_world-adult-healthy.png", dpi=300, height=6, width=12)

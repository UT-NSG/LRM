# In the Name of God
# This script allows the user to plot a multiplot entity from multiple aggregated assessment data structure 


# install.packages("ggpubr")
# setwd("~/GitHub/utlp-spm-mux")
library(ggplot2) 
library(ggpubr)
library(igraph)
library(LinkPrediction)
library(lattice)

#### Loading Saved Files ####

assessment_1 = 
  readRDS(file = "RDS_Data/AggregatedAssesment_Air-Train-Train-k=1-paper.rds")

assessment_2 = 
  readRDS(file = "RDS_Data/AggregatedAssesment_Air-Train-Train-k=69-paper.rds")

assessment_3 = 
  readRDS(file = "RDS_Data/AggregatedAssesment_Air-Train-Train-k=69-paper.rds")

AggregatedAssesment = 
  bind_rows(
    assessment_1,
    assessment_2,
    assessment_3
  )
#### Creating Plots Using Link prediction Function ####
# Metric is determined here
plot_lp_performance(AggregatedAssesment, metric = "auroc", err = "sd")
# metrics include "auroc" "recall_at_k" "avg_prec"

saveRDS(AggregatedAssesment, file = 
          "RDS_Data/MergedAssesment_Air-Train-Air-k=30-MLRM.rds"
)
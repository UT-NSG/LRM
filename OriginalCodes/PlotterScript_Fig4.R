# In the Name of God
# This script allows the user to plot a multiplot entity from multiple aggregated assessment data structure 


# install.packages("ggpubr")
#setwd("C:/Users/Amir Hekmat/Dropbox/Projects/Link Prediction/utlp-spm-mux/Old Codes")
setwd("C:/Users/Fava/Dropbox/utlp-spm-mux/Old Codes")
library(ggplot2) 
library(ggpubr)
library(igraph)
library(LinkPrediction)
library(lattice)

#To support new graphing capabilities 
source('RevisedPlotter.R')

#### Loading Saved Files ####

assessment_11 = 
  readRDS(file = "Fig4_Revised_paper.rds")

assessment_12 = 
  readRDS(file = "Fig4_Revised_paper.rds")


#### Creating Plots Using Link prediction Function ####
# Metric is determined here
# Version could be "FracRem" or "NumEig"
P11 = RevisedPlotter(assessment_11, metric = "auroc", err = "sd", Version = "NumEig")
P12 = RevisedPlotter(assessment_12, metric = "recall_at_k", err = "sd", Version = "NumEig")

# Changing labels here is possible
# y changing y labels 
# To add k=1 \n is used to go to the next line
P11$labels$y = "AUROC"
P12$labels$y = "Precision"

# omiting x label
# P11$labels$x = ""
# P12$labels$x = ""
# P13$labels$x = ""
# P21$labels$x = ""
# P22$labels$x = ""
# P23$labels$x = ""

# Adding title to the upper row
# theme(plot.title = element_text(hjust = 0.5) places the title in the middel
P11 = P11 + ggtitle("Brain/Structure \n Fraction of removed links = 0.1") +  theme(plot.title = element_text(hjust = 0.5))
P12 = P12 + ggtitle("Brain/Structure \n Fraction of removed links = 0.1") +  theme(plot.title = element_text(hjust = 0.5))

# deactivating background grid 
#P11$theme$panel.grid = element_blank()
#P12$theme$panel.grid = element_blank()

# changing ticks on x axix. 
# seq(0.0, 1, 0.2) means start from 0 to 1 by 0.2 steps.
# limits defines start and end.
# P11 = P11  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
# P12 = P12  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))

# changing ticks on y axix. 
# seq(0.0, 1, 0.2) means start from 0 to 1 by 0.2 steps.
# limits defines start and end.
# to change start and end, both breaks=seq(0,1,by=0.2), limits=c(0,1) should be changed
# for example breaks=seq(0.5,1,by=0.1), limits=c(0.5,1) changes limits from 0.5 to 1 by 0.1 steps 
P11 = 
  P11 + 
  scale_y_continuous(breaks=seq(0.5, 1, by=0.1), 
                     limits=c(0.5, 1)) 
P12 = 
  P12 + 
  scale_y_continuous(breaks=seq(0, 0.5, by=0.1), 
                     limits=c(0, 0.5)) 

# emits the legend title which used to be "predictor"
P11$theme$legend.title = element_blank()
P12$theme$legend.title = element_blank()

# It legend are to be omitted, the following lines should be uncommented
# P11$theme$legend.position = 'none'
# P12$theme$legend.position = 'none'


# drawing the multiline 
# ncol determines the number of columns
# nrow determines the number of rows
ggarrange(P11, 
          ncol = 1,
          nrow = 1)

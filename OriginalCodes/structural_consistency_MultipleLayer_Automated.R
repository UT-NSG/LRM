# In the Name of God

# This scripts show the case in which the ComparingMethodsCallerFunction is used

#### Preparations ####

# Clearing the workspace 
remove(list = ls())
# Start timing
start_time = Sys.time()
# Setting working directory to access functions and data sets
setwd("~/GitHub/utlp-spm-mux")

# Adding needed libraries 
library(igraph)
library(LinkPrediction)
library(lattice)

# Adding needed source codes 
# Creates adjacency matrix from edge list 
source("form_net.R")
# assesment function that is capable of handling other layers as input
source("prune_recover_MultipleLayer.R")

source("lp_spm_mux_coeffs.R")
source('structural_similarity.R')
source('lp_rnd.R')
source('lp_mux_spm.R')
source('lp_spm_custom.R')

source('LinkPrediction_Algorithm0.R')

source('ComparingMethodsCallerFunction.R')

source('structural_consistency_MultipleLayer.R')

#### Preparing Data Sets ####
network_layer1_name <- 'hb_l1'
network_layer2_name <- 'hb_l2'
# network_layer3_name <- 'LondonTransport_l3'

net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
# net_l3_fpath <- paste0('mux-data/', network_layer3_name, '.txt')


# network_layer1_name <- 'aarhus_l1'
# network_layer2_name <- 'aarhus_l2'
# network_layer3_name <- 'aarhus_l3'
# network_layer4_name <- 'aarhus_l4'
# network_layer5_name <- 'aarhus_l5'
# 
# net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
# net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
# net_l3_fpath <- paste0('mux-data/', network_layer3_name, '.txt')
# net_l4_fpath <- paste0('mux-data/', network_layer4_name, '.txt')
# net_l5_fpath <- paste0('mux-data/', network_layer5_name, '.txt')

# The main layer to be predicted is always layer A. 
# All other layers act as auxiliary 
gA <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)
gB <- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
# gC <- graph_from_adjacency_matrix(form_net(net_l3_fpath), mode = "undirected", diag = FALSE)
# gD <- graph_from_adjacency_matrix(form_net(net_l4_fpath), mode = "undirected", diag = FALSE)
# gE <- graph_from_adjacency_matrix(form_net(net_l5_fpath), mode = "undirected", diag = FALSE)
# Randomized graphes
gRA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gA, type = "both")), graph.density(gA), type = c("gnp"), directed = FALSE, loops = FALSE)
gRB <<-erdos.renyi.game(nrow(as_adjacency_matrix(gB, type = "both")), graph.density(gB), type = c("gnp"), directed = FALSE, loops = FALSE)
# gRC <<-erdos.renyi.game(nrow(as_adjacency_matrix(gC, type = "both")), graph.density(gC), type = c("gnp"), directed = FALSE, loops = FALSE)

#### Analysis ####

Target_Layer_Graph = gB
ReconstructionList_Graph = 
  list( 
    # list(gB), 
    # list(gC), 
    # list(gD), 
    # list(gRC),
    list(gA)
  )

structural_consistency_MultipleLayer(gB, k = 10, p_H = 0.9, auxiliaryLayers = ReconstructionList_Graph)
mean(structural_consistency(gB), k = 10, p_H = 0.9)
end_time = Sys.time()
end_time - start_time
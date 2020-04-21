library(igraph)
library(LinkPrediction)
library(lattice)

#### Preparations ####
start_time <- Sys.time()
## Set the working directory to access source files and data.
# setwd('/home/hossein/source/utlp/spm-mux/')
# setwd('C:/Users/Amir Hekmat/Documents/GitHub/utlp-spm-mux')
# setwd('C:/Users/Fava/Desktop/Dars/spm-mux/')

source('structural_similarity.R')
# 
#### Adding source files ####
source('lp_rnd.R')
source('lp_mux_spm.R')

## SPM + Our Algorithm 0
source('lp_spm_mux_coeffs.R')

## SPM + Our Algorithm 0, with correction of A changed to AoM
source('lp_spm_mux_coeffs_A_changed_to_AoM.R')

## SPM method, with custom number of eigenvalues
source('lp_spm_custom.R')

## Our Algorithm 1.2. This algorithm receives should receive AoM, M, and B
source('lp_B_Algorithm_1_2.R')

## Our Algorithm 1.1. This algorithm receives should receive AoM, M, and a mute ratio.
## This function uses A, which is not sampled, to calculate the A_tilde
source('lp_Algorithm_1_1_Using_A.R')

## Our Algorithm 1.1. This algorithm receives should receive AoM, M, and a mute ratio.
## This function uses A_R, which is sampled, to calculate the A_tilde
source('lp_Algorithm_1_1_Using_A_R.R')

source("form_net.R")

source("prune_recover_custom.R")
source("prune_recover_with_nolinks.R")
source("lp_spm_with_nolinks.R")

#### Reading Data ####

## Define Names of Layers in the Dataset
## at: air-train dataset, connection between cities by airplane and train.
## layer names: 
## physicians_l1+physicians_l2, 
## physicians_l2+physicians_l3
## physicians_l1+physicians_l3 
## at_l1+at_l2, 
## celegans_l1+celegans_l2
## celegans_l2+celegans_l3,
## celegans_l1+celegans_l3
## arxiv_l1+arxiv_l2
## arxiv_l4+arxiv_l2
## arxiv_l4+arxiv_l1
## arxiv_l2+arxiv_l8
## arxiv_l2+arxiv_l6
## arxiv_l5+arxiv_l2
## arxiv_l1+arxiv_l5
## arxiv_l2+arxiv_l7
## arxiv_l3+arxiv_l4
## hb_l1+hb_l2
## sp_l1+sp_l2
## sp_l1+sp_l3
## sp_l1+sp_l4
## sp_l1+sp_l5
## sp_l3+sp_l5
## sp_l2+sp_l4
## sp_l2+sp_l3
## sp_l3+sp_l4
## rat_l1+rat_l2
## dm_l1+dm_l2

network_layer1_name <- 'dm_l1'
network_layer2_name <- 'dm_l2'

## complete the name of path, in which the dataset is stored. 
net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')


#### Initializing the global parameters ####

## Reading the text file of data, and turning it into a graph. 
## This graph represents Matrix A in our notation.
# The graph form is a list, that shows connection between different nodes. 
# graphT_table <- read.table(paste0('mux-data/', network_layer1_name, '.txt'))
# gT <- graph.data.frame(graphT_table[, c(1, 2)], directed = FALSE)
# gT <- simplify(gT, remove.multiple = TRUE, remove.loops = TRUE)
gT <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)

## Reading the text file of data, and turning it into a graph
## This graph represents Matrix B in our notation.
# graphA_table <- read.table(paste0('mux-data/', network_layer2_name, '.txt'))
# gA <- graph.data.frame(graphA_table[, c(1, 2)], directed = FALSE)
# global_gA <<- simplify(gA, remove.multiple = TRUE, remove.loops = TRUE)
global_gA <<- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
# g <- jazz_collab

# ADJ <- as_adjacency_matrix(jazz_collab, type = "both")
# XLX <- eigen(ADJ, symmetric = T)

# Calc structural consistency
# sigma_c <- structural_consistency(gT)
# print(paste('Structural Consistency of', network_layer1_name, ':', mean(sigma_c)))

# source('spm_correlation.R')
# # spm_correlation(gT)
# mean(replicate(100, spm_correlation(gT)))
# 
# # 

## Percentages of permutated indices, in order to create A_perturbed
global_mux_pH <<- 0.1
# source('spm_correlation_mux.R')
# mean(replicate(100, spm_correlation_mux(gT)))

## Number of biggest eigenvalues and their corresponsding eigenvectors to reconstruct A_tilde
global_k <<- 10

# gT <<-erdos.renyi.game(69, 0.3, type = c("gnp"), directed = FALSE, loops = FALSE)
# global_gA <<-erdos.renyi.game(69, 0.5, type = c("gnp"), directed = FALSE, loops = FALSE)
SS <- structural_similarity(gT)
levelplot(t(SS))

# assessment <- prune_recover(gT, "lp_rnd", "lp_spm", "lp_mux_spm", "lp_spm_custom")
# assessment <- prune_recover(gT, "lp_mux_spm", "lp_spm_custom")

# assessment <- prune_recover_custom(gT, "lp_spm_mux_coeffs", "lp_spm_custom")

#### Evaluating different methods ####
## This function assess and compare the two methods given to it
# assessment <- prune_recover_custom(gT,
#               "lp_spm_mux_coeffs",
#               "lp_spm_custom",
#               "lp_B_Algorithm_1_2")

# assessment <- prune_recover_custom(gT,
#             "lp_spm_custom",
#             "lp_spm_mux_coeffs",
#             "lp_Algorithm_1_1_Using_A",
#             "lp_Algorithm_1_1_Using_A_R",
#             "lp_B_Algorithm_1_2")

# assessment <- prune_recover_custom(gT, "lp_spm_mux_coeffs","lp_B_Algorithm_1_2","lp_spm_custom","lp_Algorithm_1_1_Using_A","lp_Algorithm_1_1_Using_A_R")

assessment <- prune_recover_with_nolinks(gT, "lp_spm_with_nolinks")



# Ploting the results -------------------------------------------------------------
# metric: recall_at_k, aupr, auroc, avg_prec
plot_Handle <-plot_lp_performance(assessment, metric = "auroc", err = "sd")
# plot_Handle <- plot_lp_performance(assessment, metric = "auroc", err = "sd")
# plot_Handle$data$method[1:9]="SPM"
# plot_Handle$data$method[10:18]="Alg 1.1 A"
# plot_Handle$data$method[19,27]="Alg 1.1 A_R"
# plot_Handle
# plot_Handle$data$method[1:9]="Alg1.1A"
# plot_Handle$data$method[10:18]="Alg1.1AR"
# plot_Handle$data$method[19:27]="Alg1.2B"
# plot_Handle$data$method[28:36]="SPM"
# plot_Handle$data$method[37:45]="SPM+Alg0"
plot_Handle
end_time <- Sys.time()
print(end_time - start_time)



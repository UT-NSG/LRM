remove(list = ls())
setwd("~/GitHub/utlp-spm-mux")
# setwd("C:/Users/Fava/Desktop/Github/utlp-spm-mux")

library(igraph)
library(LinkPrediction)
library(lattice)
source("form_net.R")
source("prune_recover_custom.R")
source("lp_spm_mux_coeffs.R")
source('structural_similarity.R')
source('structural_similarity_measure_function.R')
source('structural_similarity_measure_function_Rewised_NoPlot.R')

NumberOfIterationsToAverageOver = 50
PercentOfLeading = 1

network_list_1 <- c(
  'physicians_l1',
  'physicians_l1',
  'physicians_l2',
  'aarhus_l1',
  'aarhus_l1',
  'aarhus_l1',
  'aarhus_l1',
  'aarhus_l2',
  'aarhus_l2',
  'aarhus_l2',
  'aarhus_l3',
  'aarhus_l3',
  'aarhus_l4',
  'hb_l1',
  'celegans_l1',
  'celegans_l1',
  'celegans_l2',
  'dm_l1',
  'at_l1',
  'LondonTransport_l1',
  'LondonTransport_l1',
  'LondonTransport_l2'
  )
network_list_2 <- c(
  'physicians_l2',
  'physicians_l3',
  'physicians_l3',
  'aarhus_l2',
  'aarhus_l3',
  'aarhus_l4',
  'aarhus_l5',
  'aarhus_l3',
  'aarhus_l4',
  'aarhus_l5',
  'aarhus_l4',
  'aarhus_l5',
  'aarhus_l5',
  'hb_l2',
  'celegans_l2',
  'celegans_l3',
  'celegans_l3',
  'dm_l2',
  'at_l2',
  'LondonTransport_l2',
  'LondonTransport_l3',
  'LondonTransport_l3'
  )

NumberOfPairs = 22

# Measure is L-L, L-R, R-L, R-R, pval-L-R, pval-R-L, pval-R-R
# All random substitutes have the same density as the layer 
Measure_Matrix = matrix(0, nrow = NumberOfPairs, ncol = 7)

for (LayerPair in 1:NumberOfPairs) {
  print(LayerPair)
  network_layer1_name <- network_list_1[LayerPair]
  network_layer2_name <- network_list_2[LayerPair]
  
  net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
  net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
  
  gT <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)
  global_gA <<- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
  
  # creating a same size graph by adding empty vertices
  if(vcount(gT) > vcount(global_gA)) {
    global_gA = add_vertices(global_gA, vcount(gT) - vcount(global_gA))
  }
  if(vcount(gT) < vcount(global_gA) ) {
    gT = add_vertices(gT, vcount(global_gA) - vcount(gT))
  }
  
  NCGT = nrow(as_adjacency_matrix(gT, type = "both"))
  
  # print(NCGT)
  
  DENSITY_GT = sum(as_adjacency_matrix(gT, type = "both"))/(NCGT*NCGT-NCGT)
  
  # print(DENSITY_GT)
  
  NCGA = nrow(as_adjacency_matrix(global_gA, type = "both"))
  
  # print(NCGA)
  
  DENSITY_GA = sum(as_adjacency_matrix(global_gA, type = "both"))/(NCGA*NCGA-NCGA)
  
  # print(DENSITY_GA)
  
  # gT <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GT, type = c("gnp"), directed = FALSE, loops = FALSE)
  # global_gA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GA, type = c("gnp"), directed = FALSE, loops = FALSE)
  
  # L-L
  SS <- structural_similarity(gT)
  # Me <-structural_similarity_measure_function(SS)
  Me <- structural_similarity_measure_function_Rewised_NoPlot(SS, PercentOfLeading)
  Measure_Matrix[LayerPair,1]=Me
  
  # L-R
  AcMe = 0
  AcSM = 0
  for (iir in 1:NumberOfIterationsToAverageOver){
    global_gA <<-erdos.renyi.game(nrow(as_adjacency_matrix(global_gA, type = "both")), DENSITY_GA, type = c("gnp"), directed = FALSE, loops = FALSE)
    # creating a same size graph by adding empty vertices
    if(vcount(gT) > vcount(global_gA)) {
      global_gA = add_vertices(global_gA, vcount(gT) - vcount(global_gA))
    }
    if(vcount(gT) < vcount(global_gA) ) {
      gT = add_vertices(gT, vcount(global_gA) - vcount(gT))
    }
    SS <- structural_similarity(gT)
    # Me <-structural_similarity_measure_function(SS)
    Me <- structural_similarity_measure_function_Rewised_NoPlot(SS, PercentOfLeading)
    AcMe = AcMe + Me
    AcSM = AcSM + Me*Me
  }
  AvR = AcMe/NumberOfIterationsToAverageOver
  Measure_Matrix[LayerPair,2]= AvR
  Var = (AcSM )/NumberOfIterationsToAverageOver - AvR*AvR
  SDD = sqrt(Var) 
  Measure_Matrix[LayerPair,3]= 1 - pnorm((Measure_Matrix[LayerPair,1] - AvR)/SDD)
  # R-L
  AcMe = 0
  AcSM = 0
  for (iir in 1:NumberOfIterationsToAverageOver){
    global_gA <<- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
    gT <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GT, type = c("gnp"), directed = FALSE, loops = FALSE)
    # creating a same size graph by adding empty vertices
    if(vcount(gT) > vcount(global_gA)) {
      global_gA = add_vertices(global_gA, vcount(gT) - vcount(global_gA))
    }
    if(vcount(gT) < vcount(global_gA) ) {
      gT = add_vertices(gT, vcount(global_gA) - vcount(gT))
    }
    SS <- structural_similarity(gT)
    # Me <-structural_similarity_measure_function(SS)
    Me <- structural_similarity_measure_function_Rewised_NoPlot(SS, PercentOfLeading)
    AcMe = AcMe + Me
    AcSM = AcSM + Me*Me
  }
  AvR = AcMe/NumberOfIterationsToAverageOver
  Measure_Matrix[LayerPair,4]= AvR
  Var = (AcSM )/NumberOfIterationsToAverageOver - AvR*AvR
  SDD = sqrt(Var) 
  Measure_Matrix[LayerPair,5]= 1 - pnorm((Measure_Matrix[LayerPair,1] - AvR)/SDD)
  
  # R-R
  AcMe = 0
  AcSM = 0
  for (iir in 1:NumberOfIterationsToAverageOver){
    gT <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GT, type = c("gnp"), directed = FALSE, loops = FALSE)
    global_gA <<-erdos.renyi.game(nrow(as_adjacency_matrix(global_gA, type = "both")), DENSITY_GA, type = c("gnp"), directed = FALSE, loops = FALSE)
    # creating a same size graph by adding empty vertices
    if(vcount(gT) > vcount(global_gA)) {
      global_gA = add_vertices(global_gA, vcount(gT) - vcount(global_gA))
    }
    if(vcount(gT) < vcount(global_gA) ) {
      gT = add_vertices(gT, vcount(global_gA) - vcount(gT))
    }
    SS <- structural_similarity(gT)
    # Me <-structural_similarity_measure_function(SS)
    Me <- structural_similarity_measure_function_Rewised_NoPlot(SS, PercentOfLeading)
    AcMe = AcMe + Me
    AcSM = AcSM + Me*Me
  }
  AvR = AcMe/NumberOfIterationsToAverageOver
  Measure_Matrix[LayerPair,6]= AvR
  Var = (AcSM )/NumberOfIterationsToAverageOver - AvR*AvR
  SDD = sqrt(Var) 
  Measure_Matrix[LayerPair,7]= 1 - pnorm((Measure_Matrix[LayerPair,1] - AvR)/SDD)
}


print(Measure_Matrix)

RName = paste(network_list_1,network_list_2)
DF_Measure = data.frame(Measure_Matrix, row.names = RName)  
colnames(DF_Measure) <- c("LL", "LR", "pValue LR", "RL", "pValue RL", "RR", "pValue RR")
write.table(DF_Measure, "MeasureMatrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)




# In the Name of God
# Comparing Methods R Function
# Question: Why use 'LA' 
# in calculation of first k eigenvalues
# 'LA' mean largest algebraically
# is there a reason? //???????

#### Description ####
# Changed all variables to local. No global is left. 
# Also added comments to make the code more readable

# Version Dates:
# V1: 13980224
# V2: 13980225
# V3: 13980301

#### Preparations ####
# Adding needed libraries 
library(tibble)

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
# It is advised to pass same size graphes to the function
ComparingMethodsCallerFunctionWrapper <- 
  function (Target_Layer_Graph, 
            method = "LinkPrediction_Algorithm0", 
            ReconstructionList_Graph = NA, 
            LegendList,
            probes = seq(0.1, 0.9, 0.1), 
            epochs = 30,
            PowerSelector = NA,
            NumberOfLeadingEigenValuesVec = c(1,2,3,5, 10, 20),
            SelectorOfEigenValues = "LA") 
    {
    AggregatedAssesmentAggregated = NA
    for(iC in 1:length(NumberOfLeadingEigenValuesVec)) {
      NumberOfLeadingEigenValues = NumberOfLeadingEigenValuesVec[iC]
      print(NumberOfLeadingEigenValues)
      AggregatedAssesment = 
        ComparingMethodsCallerFunction(Target_Layer_Graph,
        method, 
        ReconstructionList_Graph, 
        LegendList,
        probes, 
        epochs,
        PowerSelector,
        NumberOfLeadingEigenValues,
        SelectorOfEigenValues)
      
      AggregatedAssesment = add_column(AggregatedAssesment, k = NumberOfLeadingEigenValues)
      
      if(is.na(AggregatedAssesmentAggregated)) {
        AggregatedAssesmentAggregated = AggregatedAssesment;
      }
      else{
        AggregatedAssesmentAggregated = bind_rows(AggregatedAssesmentAggregated, AggregatedAssesment)
      }
    
      
      }
    
    
    
    return(AggregatedAssesmentAggregated)
}


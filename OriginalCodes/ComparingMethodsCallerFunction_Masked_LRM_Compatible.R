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
library(igraph)
library(LinkPrediction)
library(lattice)

# Adding needed source codes 
# Creates adjacency matrix from edge list 
source("form_net.R")
# assesment function that is capable of handling other layers as input
source("prune_recover_MultipleLayer_Masked_LRM_Compatible.R")
source("lp_spm_mux_coeffs.R")
source('structural_similarity.R')
source('lp_rnd.R')
source('lp_mux_spm.R')
source('lp_spm_custom.R')

# It is advised to pass same size graphes to the function
ComparingMethodsCallerFunction_Masked_LRM_Compatible <- 
  function (Target_Layer_Graph, 
            method, 
            ReconstructionList_Graph = NA, 
            LegendList,
            includeSelf = 1,
            probesOfLinks = seq(0.1, 0.9, 0.1), 
            NoLinksToLinksRatio=2,
            epochs = 30,
            PowerSelector = NA,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA") 
    {
    
    # Masked LRM using only self Method
    if(includeSelf == 1){
      assessment_MLRM <- 
        prune_recover_MultipleLayer_Masked_LRM_Compatible(
          Target_Layer_Graph,
          method,
          probesOfLinks = probesOfLinks, 
          NoLinksToLinksRatio = NoLinksToLinksRatio,
          epochs = epochs, 
          preserve_conn = FALSE, 
          use_weights = FALSE, 
          includeSelf = 1,
          auxiliaryLayers = NA, # NA means SPM, auxiliaryLayers uses these layers
          PowerSelector = PowerSelector,
          NumberOfLeadingEigenValues = NumberOfLeadingEigenValues,
          SelectorOfEigenValues = SelectorOfEigenValues
        )
      # Changing label
      assessment_MLRM$method = "MLRM(Self)" #?????CP
      AggregatedAssesment = assessment_MLRM
    }
    else
    {
      AggregatedAssesment = NA;
    }
    
    if(is.na(ReconstructionList_Graph)==FALSE)
    {
      # For each Graph List
      for(i_c_g in 1:length(ReconstructionList_Graph)) {
        auxiliaryLayers = ReconstructionList_Graph[[i_c_g]]
        NumAuxLay = length(auxiliaryLayers)
        
        # Making sure that all graphes are the same size 
        maximumVcount = vcount(Target_Layer_Graph)
        for (i_m_c in 1:NumAuxLay) {
          maximumVcount = max(maximumVcount, vcount( auxiliaryLayers[[i_m_c]] ) )
        }
        if(vcount(Target_Layer_Graph) < maximumVcount) {
          Target_Layer_Graph = 
            add_vertices(Target_Layer_Graph, maximumVcount - vcount(Target_Layer_Graph))
        }
        for (i_m_v in 1:NumAuxLay) {
          if(vcount( auxiliaryLayers[[i_m_v]]) < maximumVcount) {
            auxiliaryLayers[[i_m_v]] = 
              add_vertices(auxiliaryLayers[[i_m_v]], 
                           maximumVcount - vcount(auxiliaryLayers[[i_m_v]]))
          }
        }
        #calling the assessment function for each set 
        assessment_AUX <- 
          prune_recover_MultipleLayer_Masked_LRM_Compatible(
            Target_Layer_Graph,
            method,
            probesOfLinks = probesOfLinks, 
            NoLinksToLinksRatio = NoLinksToLinksRatio,
            epochs = epochs, 
            preserve_conn = FALSE, 
            use_weights = FALSE, 
            includeSelf = includeSelf,
            auxiliaryLayers = auxiliaryLayers, # NA means SPM, auxiliaryLayers uses these layers
            PowerSelector = PowerSelector,
            NumberOfLeadingEigenValues = NumberOfLeadingEigenValues,
            SelectorOfEigenValues = SelectorOfEigenValues
          )
        assessment_AUX$method = LegendList[[i_c_g]]
        #appending the result to the previous ones 
        if(is.na(AggregatedAssesment)==TRUE)
        {
          AggregatedAssesment = assessment_AUX;
        }
        else
        {
          AggregatedAssesment = 
            bind_rows(AggregatedAssesment,
                      assessment_AUX)
        }
        
      }
    }
    
  
    return(AggregatedAssesment)
}


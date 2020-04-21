# In the Name of God

# This scripts show the case in which the ComparingMethodsCallerFunction is used

#### Preparations ####

# Clearing the workspace 
remove(list = ls())
# Start timing
start_time = Sys.time()
# Setting working directory to access functions and data sets
# setwd("C:/Users/dozha/Dropbox/utlp-spm-mux/Old Codes")
# setwd("C:/Users/Fava/Dropbox/utlp-spm-mux/Old Codes")
setwd("C:/Users/Amir Hekmat/Dropbox/Projects/Link Prediction/utlp-spm-mux/Old Codes")

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
source('LinkPrediction_Algorithm0_OnlyLRMSupport.R')
source('LinkPrediction_Algorithm0_Merge_By_GOR.R')
source('LinkPrediction_Algorithm0_Merge_By_PCC.R')
source('LinkPrediction_Algorithm0_Merge_By_Norm_GOR.R')
source('LinkPrediction_Algorithm0_Merge_By_Norm_PCC.R')



source('ComparingMethodsCallerFunction.R')
source('ComparingMethodsCallerFunctionWrapper.R')
#To support new graphing capabilities 
source('RevisedPlotter.R')

NumberOfLeadingEigenValues_Vector = c(1, 5, 10, 15, 20, 25)

Duplex_Table = read.table('data_spec_DrivedFromGoogleSheet.txt', sep = ',', stringsAsFactors=FALSE)
nM = length(NumberOfLeadingEigenValues_Vector)
Results_Matrix = matrix(
  0, 
  ncol = nM, 
  nrow = nrow(Duplex_Table) 
)

Results_Matrix_SD = matrix(
  0, 
  ncol = nM, 
  nrow = nrow(Duplex_Table) 
)

Name_Of_Rows = c()

for (iNR in 1:nrow(Duplex_Table) ) {
  
  network_layer_name_Tar = Duplex_Table[iNR, 1]
  network_layer_name_Aux = Duplex_Table[iNR, 2]
  
  # forming the path string
  net_l_fpath_Tar = 
    paste0('../mux-data/', network_layer_name_Tar, '.txt')
  net_l_fpath_Aux = 
    paste0('../mux-data/', network_layer_name_Aux, '.txt')
  # loading data from file and forming graph objects
  Target_Layer_Graph = 
    igraph::graph_from_adjacency_matrix(
      form_net(net_l_fpath_Tar), 
      mode = "undirected", diag = FALSE)
  gAux = 
    igraph::graph_from_adjacency_matrix(
      form_net(net_l_fpath_Aux), 
      mode = "undirected", diag = FALSE)
  
  Name_Of_Row = paste0(network_layer_name_Tar,',',network_layer_name_Aux)
  Name_Of_Rows = rbind(Name_Of_Rows, Name_Of_Row)
  
  for (iEV in 1:nM) {
    NoLEV = NumberOfLeadingEigenValues_Vector[iEV]
    
    print(paste0('Analyzing Row ', Name_Of_Row, ', k=', NoLEV))
    ReconstructionList_Graph = list(list(gAux))
    LegendList =list("SPM(Self)+GOR+LRM(AUX)")
    
    Assessed =  tryCatch({
      Assessed =
        ComparingMethodsCallerFunctionWrapper(
          Target_Layer_Graph, 
          method = "LinkPrediction_Algorithm0_Merge_By_Norm_PCC", 
          ReconstructionList_Graph = ReconstructionList_Graph, 
          LegendList,
          probes = seq(from = 0.1, to = 0.1, by = 0.1), 
          epochs = 30, 
          PowerSelector = NA,
          NumberOfLeadingEigenValuesVec = c(NoLEV), 
          SelectorOfEigenValues = "LA"
        )
    }, 
    error = function(err) {
      # error handler picks up where error was generated
      print(paste("An error occured:  ",err))
      print('Using zero as score')
      Assessed$auroc = 0
      
      return(Assessed)
    }, finally = {
      
    }) # END tryCatch
    Results_Matrix[iNR, iEV] = mean(Assessed$auroc)
    Results_Matrix_SD[iNR, iEV] = sd(Assessed$auroc)
  }
}

dimnames(Results_Matrix) = list(Name_Of_Rows)
dimnames(Results_Matrix_SD) = list(Name_Of_Rows)
write.table(Results_Matrix, 'SPM(Self)LRM(AUX)NormGOR_mean.txt', row.names = TRUE)
# proper naming scheme
# atl1_atl2_13981008_mean.txt
write.table(Results_Matrix_SD, 'SPM(Self)LRM(AUX)NormGOR_std.txt', row.names = TRUE)
# proper naming scheme
# atl1_atl2_13981008_std.txt
print(Results_Matrix)
print(Results_Matrix_SD)
# Showing the elapsed time
end_time = Sys.time()
end_time - start_time




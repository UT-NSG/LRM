# In the Name of God

## This files calls the methods defined in "AnalyzeKey" 
## on the datasets defined in the file "data_spec_handson.txt"
## and writes the answeres to the file "ResutlsTextFile.txt".


# This scripts calculates the AUC automatically for a number of multiplex networks, defined by data_spec.txt file
# Right now, only duplex networks are supprted. ###To all later, multiple Auxilary layers###
#### Preparations ####
# Clearing the workspace 
remove(list = ls())

# Start timing
start_time = Sys.time()

# Setting working directory to access functions and data sets
setwd("C:/Users/Amir Hekmat/Dropbox/Projects/Link Prediction/utlp-spm-mux/LP New Code")
# setwd("C:/Users/Fava/Dropbox/utlp-spm-mux/LP New Code")

## Adding needed libraries 
# to use graph_from_adjacency_matrix and erdos.renyi graph
library(igraph)

# to use plot_lp plotting function
library(LinkPrediction)


# library(lattice)
# library(matrixcalc) # install.packages("matrixcalc")

## Adding needed source codes 
# Creates adjacency matrix from edge list 
source('Core/form_net.R')

# Main function that assesses the method

source('Core/N_Assess_Method.R')

# List of methods to be used on the ListOfSources sources. It should be the same length.
# New methods are all preceded by NS_
# current options: 'NS_SPM', 'NS_LRM', 'NS_LRM_Perturbed', 'NS_Masked_LRM', 'NS_Masked_LRM_SmartNumEigen'

NumberOfLeadingEigenValues_Vector = 
  c(1, 5, 10, 15, 20, 25, 30)
  # c(1, 5, 10)
  # c(29, 30, NA)
  # c(1, 5)
  # seq(1, 10, 3)


# Method by which the scores of different methods are aggregated
# current options: 'SimpleAddition'
# 'SimpleAddition' adds the score of different items in the list together
# Supported Methods Include: 'SimpleAddition', 'PCC', 'GOR'
# AggregationMethod = 'GOR' # the lost jewel of link prediction
# Now the aggregation method is also part of the AnalyzeKey


# Extra processing done on the data
# current options: 'None', 'Laplacian', 'Laplacian_A', 'Normalized_Laplacian', 'Normalized_Laplacian_A'
# 'None' does no extra processing
# 'Laplacian' uses laplacian matrix of the graph for reconstruction, based on adjacency reconstruction 
# 'Laplacian_A' uses laplacian matrix of the graph for reconstruction, based on laplacian reconstruction
# 'Normalized_Laplacian' normalized version of Laplacian
# 'Normalized_Laplacian_A' normalized version of Laplacian_A

ExtraProcessing = 'None'



# AUC of the following modes will be calculated for each duplex
# Any combination following the model is possible.
AnalyzeKey =
  list(
    list(
      list('Train'),
      list('NS_SPM'),
      'SimpleAddition',
      0
    )
  )
# AnalyzeKey =
#   list(
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_LRM'),
#       'SimpleAddition',
#       0
#     ),
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_LRM_Perturbed'),
#       'SimpleAddition',
#       0
#     ),
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_LRM_Perturbed'),
#       'GOR',
#       0
#     ),
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_SPM'),
#       'SimpleAddition',
#       0
#     )
#   )
# AnalyzeKey =
#   list(
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_LRM_Perturbed'),
#       'GOR',
#       0
#     ),
#     list(
#       list('Train',  'AUX'),
#       list('NS_SPM', 'NS_LRM_Perturbed'),
#       'PCC',
#       0
#     )
#   )
# AnalyzeKey = 
#   list(
#     list(
#       list('Train',  'AUX'), 
#       list('NS_SPM', 'NS_SPM')
#     ),
#     list(
#       list('Train',  'AUX'), 
#       list('NS_SPM', 'NS_Masked_LRM')
#     ),
#     list(
#       list('Train',  'AUX'), 
#       list('NS_SPM', 'NS_Masked_LRM_SmartNumEigen')
#     )
#   )


Duplex_Table = read.table('data_spec.txt', sep = ',', stringsAsFactors=FALSE)

nM = length(NumberOfLeadingEigenValues_Vector)

Results_Matrix = matrix(
  0, 
  ncol = length(AnalyzeKey), 
  nrow = nM*nrow(Duplex_Table) 
  )

Results_Matrix_SD = matrix(
  0, 
  ncol = length(AnalyzeKey), 
  nrow = nM*nrow(Duplex_Table) 
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
  gTar = 
    igraph::graph_from_adjacency_matrix(
      form_net(net_l_fpath_Tar), 
      mode = "undirected", diag = FALSE)
  gAux = 
    igraph::graph_from_adjacency_matrix(
      form_net(net_l_fpath_Aux), 
      mode = "undirected", diag = FALSE)
  for (iEV in 1:nM) {
    NoLEV = NumberOfLeadingEigenValues_Vector[iEV]
    Name_Of_Row = paste0(network_layer_name_Tar,',',network_layer_name_Aux,', k = ',NoLEV, ': ')
    Name_Of_Rows = rbind(Name_Of_Rows, Name_Of_Row)
    
    print(paste0('Analyzing Row ', iNR, ', ', Name_Of_Row))
    
    for (iNC in 1:ncol(Results_Matrix)) { # For each key
      
      # Create list of sources 
      ListOfSources = list() 
      for (iOverK in 1:length(AnalyzeKey[[iNC]][[1]])) {
        if (AnalyzeKey[[iNC]][[1]][[iOverK]] == 'Train') {
          ListOfSources = append(ListOfSources, 'Train')
        }
        else if(AnalyzeKey[[iNC]][[1]][[iOverK]] == 'AUX') {
          ListOfSources = append(ListOfSources, list(gAux))
        }
      }
  
      # Create list of methods 
      ListOfMethods = AnalyzeKey[[iNC]][[2]]
      # print('Methods to call')
      # print(ListOfMethods)
      
      LegendStr = 'Does not matter, since there is no plot'
      
      AggregationMethod = AnalyzeKey[[iNC]][[3]]
      
      AddingAUXAdjacencyScale = AnalyzeKey[[iNC]][[4]]
      
      Assessed =  tryCatch({
        Assessed = 
          N_Assess_Method(
            Tg_Graph = gTar, 
            methodList = ListOfMethods,
            sourceList = ListOfSources,
            aggregation_method = AggregationMethod,
            extra_processing = ExtraProcessing,
            legendFig = LegendStr,
            probesOfLinks = seq(0.1, 0.1, 0.1), # Percentent of missing links  
            epochs = 30, # Number of Train and Test Sets are calculated
            NumberOfLeadingEigenValues = NoLEV, # if NA, all eigenvalues are used
            SelectorOfEigenValues = "LA", # common options include LA for algebraic largest, and LM for magnitude largest
            PowerSelector = NA, # Percent of power of eigenvalues to be used. Used by algegraic largest method. If powerselector has value, NumberOfLeadingEigenValues will be ignored. If number of eigenvalues is going to be specifically given, give NA as input for this arguement
            NoLinksToLinksRatio = 10000, # this arguement determines the ratio between no links and links
            AddingAUXAdjacencyScale = AddingAUXAdjacencyScale # This value chooses the scale of Adajency Matrix of AUX layer to be added
          )
      }, 
      error = function(err) {
        # error handler picks up where error was generated
        print(paste("An error occured big time:  ",err))
        print('Using zero as score')
        Assessed$auroc = 0

        return(Assessed)
      }, finally = {
        
      }) # END tryCatch
      Results_Matrix[(iNR-1)*nM + iEV, iNC] = mean(Assessed$auroc)
      Results_Matrix_SD[(iNR-1)*nM + iEV, iNC] = sd(Assessed$auroc)
      
    }
  }
}

dimnames(Results_Matrix) = list(Name_Of_Rows)
dimnames(Results_Matrix_SD) = list(Name_Of_Rows)
write.table(Results_Matrix, 'RawResults/SPM_AarhusL5_mean.txt', row.names = TRUE)
# proper naming scheme
# atl1_atl2_13981008_mean.txt
write.table(Results_Matrix_SD, 'RawResults/SPM_AarhusL5_std.txt', row.names = TRUE)
# proper naming scheme
# atl1_atl2_13981008_std.txt
print(Results_Matrix)
print(Results_Matrix_SD)
# Showing the elapsed time
end_time = Sys.time()
end_time - start_time

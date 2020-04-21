# for use of Matrix function
library(Matrix)

# for use of as_adjacency_matrix, graph_from_adjacency_matrix, and get_non_edges
library(igraph)

# for method caller
source('Core/NS_SPM.R')
source('Core/NS_LRM.R')
source('Core/NS_LRM_Perturbed.R')
source('Core/NS_Masked_LRM.R')
source('Core/NS_Masked_LRM_SmartNumEigen.R')

source('Core/AggregationMethodExtraParameter.R')
source('Core/SimpleAddition.R')


Extended_Method_Caller =
  function (
    g_train_Links, 
    g_train_NoLinks,
    methodList,
    sourceList,
    aggregation_method = 'SimpleAddition',
    extra_processing = 'None',
    NumberOfLeadingEigenValues = NA, 
    SelectorOfEigenValues = 'LA',
    PowerSelector = NA,
    AddingAUXAdjacencyScale = 0
  )
  {
    N = igraph::vcount(g_train_Links)
    # Make weights 1 to ensure functionality
    E(g_train_Links)$weight = 1
    if(is.na(g_train_NoLinks)==FALSE) {
      E(g_train_NoLinks)$weight = 1
    }
    
    # Observed adjacency matrix
    AoW = igraph::as_adjacency_matrix(g_train_Links, type = "both")
    W = AoW + igraph::as_adjacency_matrix(g_train_NoLinks, type = "both");
    
    # Pre extra processing
    if(extra_processing == 'None') {
      # Do extra processing here, if there is extra processing
      # preprocessing here
    }
    else if(extra_processing == 'Laplacian') {
      
    }
    else if(extra_processing == 'Laplacian_A') {
      
    }
    else if(extra_processing == 'Normalized_Laplacian') {
      
    }
    else if(extra_processing == 'Normalized_Laplacian_A') {
      
    }
    else {
      # To add other types of extra processing later
      print('Error, undefined extra processing parameter')
      print(extra_processing)
    }
    
    ScoresMatrix = Matrix::Matrix(0, N, N)
    
    # # Debugging
    # print('I am here before calling method in Extended Method Caller')
    
    for(i_s_m in 1:length(ListOfMethods)) {
      TPredMat = do.call(
        methodList[[i_s_m]],
        list(
          gLT = g_train_Links,
          gNLT = g_train_NoLinks,
          auxL = sourceList[[i_s_m]],
          NumberOfLeadingEigenValues = NumberOfLeadingEigenValues,
          SelectorOfEigenValues = SelectorOfEigenValues,
          PowerSelector = PowerSelector, 
          EPType = extra_processing
        )
      )
      
      
      if(sourceList[[i_s_m]]=='Train') {
        gAdd = g_train_Links
      }
      else {
        gAdd = sourceList[[i_s_m]]
      }
      Train_Matrix = igraph::as_adjacency_matrix(g_train_Links)
      AUX_Matrix = igraph::as_adjacency_matrix(gAdd)
      TPredMat = TPredMat + AddingAUXAdjacencyScale * AUX_Matrix
      
      # # Debuggin
      # print('I am here before calling aggregation method in Extended Method Caller')
      # print(ScoresMatrix)
      # print(TPredMat)
      
      AMEPL = AggregationMethodExtraParameter(
        aggregation_method = aggregation_method, 
        Train_Matrix, 
        AUX_Matrix, 
        EPType = extra_processing)
      
      ScoresMatrix = 
        SimpleAddition(
          ScoreSoFar = ScoresMatrix,
          NewPrediction = TPredMat,
          EPType = extra_processing,
          AggregationMethodExtraParameter = AMEPL
          )
      
      # ScoresMatrix = tryCatch({
      #   ScoresMatrix = do.call(
      #     aggregation_method,
      #     list(
      #       ScoreSoFar = ScoresMatrix, 
      #       NewPrediction = TPredMat,
      #       EPType = extra_processing,
      #       AggregationMethodExtraParameter = AMEPL
      #     )
      #   )
      # }, 
      # error = function(err) {
      #   # error handler picks up where error was generated
      #   print(paste("An error occured in the aggregation method: ",err))
      #   print('Using all zeros as score matrix')
      #   ScoresMatrix = Matrix::Matrix(0, N, N)
      #   return(ScoresMatrix)
      # }, 
      # finally = {
      #   
      # }) # END tryCatch
      
      # # Debuggin
      # print('I am here after calling aggregation method Extended Method Caller')
    }
    
    
    # Post extra processing
    if(extra_processing == 'None') {
      # Do extra processing here, if there is extra processing
      # preprocessing here
    }
    else if(extra_processing == 'Laplacian') {
      
    }
    else if(extra_processing == 'Laplacian_A') {
      # In this processing scheme, the lapalacian matrix is created, now it has to be changed into score
      ScoresMatrix =  (diag(diag(ScoresMatrix)) - ScoresMatrix)
    }
    else if(extra_processing == 'Normalized_Laplacian') {
      
    }
    else if(extra_processing == 'Normalized_Laplacian_A') {
      ScoresMatrix =  (diag(diag(ScoresMatrix)) - ScoresMatrix)
    }
    else {
      # To add other types of extra processing later
      print('Error, undefined extra processing parameter')
      print(extra_processing)
    }
    
    ################# Think about it, why does it work?
    if (SelectorOfEigenValues == 'SM' || SelectorOfEigenValues == 'SA') {
      ScoresMatrix = -ScoresMatrix
    }
    
    # # Debuggin
    # print('I am here before creating prediction')
    
    # Creating the prediction based on the ScoreMatrix. The higher the score, method predicts being link more strongly
    ObservedGraph = igraph::graph_from_adjacency_matrix(W, mode = "undirected", diag = FALSE)
    non_edges = get_non_edges(ObservedGraph)
    prediction = tibble::tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                         scr = ScoresMatrix[non_edges]) %>% dplyr::arrange(desc(scr))
    
    # # Debuggin
    # print('I am here before returning prediction in Extended Method Caller')
    
    return(prediction)
  }

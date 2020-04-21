# for E, ecount, vcount, as_adjacency_matrix, laplacian_matrix, and delete_edges
library(igraph)

# for Matrix
library(Matrix)

# for eigs_sym
library(RSpectra)


# library(dplyr)

NS_LRM_Perturbed = 
  function(
    gLT,
    gNLT,
    auxL,
    NumberOfLeadingEigenValues,
    SelectorOfEigenValues,
    PowerSelector,
    EPType = 'None',
    # Local variables to LRM_Pertured
    p_H = 0.1,
    local_epochs = 10
    )
    {
    
    if(auxL=='Train') {
      g = gLT
    }
    else {
      g = auxL
    }
    
    # Extra processing
    if(EPType == 'None') {
      # Do extra processing here, if there is extra processing
      matrix_creator_Function_Hanlge = igraph::as_adjacency_matrix
    }
    else if(EPType == 'Laplacian') {
      matrix_creator_Function_Hanlge = igraph::laplacian_matrix
    }
    else if(EPType == 'Laplacian_A') {
      matrix_creator_Function_Hanlge = igraph::laplacian_matrix
    }
    else if(EPType == 'Normalized_Laplacian') {
      matrix_creator_Function_Hanlge = function(A_Graph) {(0.5/(igraph::ecount(A_Graph)))*igraph::laplacian_matrix(A_Graph)}
    }
    else if(EPType == 'Normalized_Laplacian_A') {
      matrix_creator_Function_Hanlge = function(A_Graph) {(0.5/(igraph::ecount(A_Graph)))*igraph::laplacian_matrix(A_Graph)}
    }
    else {
      # To add other types of extra processing later
      print(EPType)
      stop('Error: Error, undefined extra processing parameter')
    }
    
    # Make weights 1 to ensure functionality
    igraph::E(g)$weight = 1
    
    # Network stats and adjacency matrix
    L = igraph::ecount(g)
    N = igraph::vcount(g)
    A = matrix_creator_Function_Hanlge(g)
    
    
    if(is.na(NumberOfLeadingEigenValues)==FALSE ) {
      k = NumberOfLeadingEigenValues
    }
    else if(is.na(NumberOfLeadingEigenValues)==TRUE) {
      k = N
    }

    AUX_Pred_Average = Matrix::Matrix(0, N, N)
    AUX_Pred_Sum = Matrix::Matrix(0, N, N)
    
    # Find the number of Eigenvalues 
    if(is.na(PowerSelector) == FALSE ){
      EIG_AUX = eigen(AUX_Matrix, symmetric = T, only.values = TRUE)
      EIG_AUX_Sorted = sort(abs(EIG_AUX$values), decreasing = TRUE)
      TotalEnergy = sum((EIG_AUX$values) * (EIG_AUX$values) );
      RequiredEnergy = TotalEnergy*PowerSelector;
      SumOfEnergySoFar = 0;
      k = 0
      while(SumOfEnergySoFar < RequiredEnergy){
        k = k+1
        SumOfEnergySoFar = SumOfEnergySoFar + EIG_AUX_Sorted[k]*EIG_AUX_Sorted[k]
      }
    }
    
    # Calculate number of links to be removed
    links_to_remove = round(L*p_H)
    
    for(i in 1:local_epochs){
      # Determine the link perturbation set dlt_E at random
      to_remove = sample(L, links_to_remove)
      gAUX_perturbed = igraph::delete_edges(g, to_remove)
      
      # making sure the size does not change
      if(igraph::vcount(gAUX_perturbed) < N) {
        gAUX_perturbed = 
          igraph::add_vertices(gAUX_perturbed, N - igraph::vcount(gAUX_perturbed))
      }
      # Compute the adjacency matrix of the graph after removal of dlt_E
      AUX_Matrix = matrix_creator_Function_Hanlge(gAUX_perturbed)

      
      AUX_Pred = Matrix::Matrix(0, N, N)
      
      
      AUX_Pred =  tryCatch({
        # Diagonalise AUX_Pred
        if(is.na(k) || k == N){
          AUX_XLX = eigen(AUX_Matrix, symmetric = T)
        }
        else{
          AUX_XLX = RSpectra::eigs_sym(AUX_Matrix, k = k, which = SelectorOfEigenValues)
        }
        
        # calculating alpha_i coefficients 
        alpha = c()
        for(i in 1:k) {
          x_i = AUX_XLX$vectors[, i]
          alpha = c(alpha, sum(diag((x_i %*% t(x_i)) %*% A)))
        }
        if(k==1){
          AUX_Pred = 
            (AUX_XLX$vectors
             %*% (alpha)
             %*% t(AUX_XLX$vectors))
        }
        else{
          AUX_Pred = 
            (AUX_XLX$vectors
             %*% diag(alpha)
             %*% t(AUX_XLX$vectors))        
        }
      }, 
      error = function(err) {
        # error handler picks up where error was generated
        print(paste("An error occured:  ",err))
        print('Using all zeros as prediction')
        AUX_Pred = Matrix::Matrix(0, N, N)
        return(AUX_Pred)
      }, finally = {
        
      }) # END tryCatch
      AUX_Pred_Sum = AUX_Pred_Sum + AUX_Pred
    }
    
    AUX_Pred_Average = AUX_Pred_Sum/local_epochs
    return(AUX_Pred_Average)
}
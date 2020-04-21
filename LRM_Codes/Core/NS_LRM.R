# for E, ecount, vcount, as_adjacency_matrix, laplacian_matrix
library(igraph)

# for Matrix
library(Matrix)

# for eigs_sym
library(RSpectra)


# library(dplyr)

NS_LRM = 
  function(
    gLT,
    gNLT,
    auxL,
    NumberOfLeadingEigenValues,
    SelectorOfEigenValues,
    PowerSelector,
    EPType = 'None'
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
    
    # for each extra layer, we calculate the score
    # If there are no extra layers, the method is just SPM
    AUX_Matrix = matrix_creator_Function_Hanlge(g)
    
    AUX_Pred = Matrix::Matrix(0, N, N)
    
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
    # print(alpha)
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

    return(AUX_Pred)
}
# for E, ecount, vcount, as_adjacency_matrix, laplacian_matrix, and delete_edges
library(igraph)

# for Matrix
library(Matrix)

# for eigs_sym
library(RSpectra)


# library(dplyr)
library(dplyr)

# install.packages('lme4')
# install.packages("matrixcalc")
# library(lme4)

library(matrixcalc)

# We should find a better name for our algorithm
NS_Masked_LRM = 
  function(gLT,
           gNLT,
           auxL,
           NumberOfLeadingEigenValues,
           SelectorOfEigenValues,
           PowerSelector,
           EPType = 'None'
           )
  {
    
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
    
    
    # compatibality of names:
    g = gLT
    CompG = gNLT
    
    # Make weights 1 to ensure functionality
    igraph::E(g)$weight = 1
    igraph::E(CompG)$weight = 1
    
    # Network stats and adjacency matrix
    L = igraph::ecount(g)
    N = igraph::vcount(g)
    AoW = matrix_creator_Function_Hanlge(g)
    W = AoW + matrix_creator_Function_Hanlge(CompG, type = "both")
    
    
    if(is.na(NumberOfLeadingEigenValues)==FALSE ) {
      k = NumberOfLeadingEigenValues
    }
    else if(is.na(NumberOfLeadingEigenValues)==TRUE) {
      k = N
    }
    
    if(auxL=='Train') {
      g = gLT
    }
    else {
      g = auxL
    }
     
    # for each extra layer, we calculate the score
    # If there are no extra layers, the method is just SPM
    AUX_Matrix = matrix_creator_Function_Hanlge(g)
    
    AUX_Pred = Matrix::Matrix(0, N, N)
    
    
    # Find the number of Eigenvalues 
    if(is.na(PowerSelector) == FALSE ){
      print("Power selector is not yet supported by NS_Masked_LRM.")
      print("Contact AHTK is implementation is required.")
    }
    
    AUX_Pred =  tryCatch({
      # Diagonalise AUX_Pred
      if(is.na(k) || k == N){
        AUX_XLX = eigen(AUX_Matrix, symmetric = T)
      }
      else{
        AUX_XLX = RSpectra::eigs_sym(AUX_Matrix, k = k, which = SelectorOfEigenValues)
      }
      
      # calculating lambda coefficients
      T_matrix = Matrix::Matrix(0, k, k);
      To_matrix = Matrix::Matrix(0, k,1);
      lambda = Matrix::Matrix(0, k, 1);
      for(i in 1:k) {
        q_i = AUX_XLX$vectors[, i]
        Q_I = q_i %*% t(q_i)
        To_matrix[i,1] = sum( diag( (Q_I*W)%*%(AoW) ) )
        for(j in 1:k) {
          q_j = AUX_XLX$vectors[, j]
          Q_J = q_j %*% t(q_j)
          T_matrix[i, j] = sum(diag( ( (Q_I)*W )%*%( (Q_J)*W )  ))
        }
      }
      
      if(is.singular.matrix(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE )
      {
        AUX_Pred = Matrix::Matrix(0, N, N)
        print('Simple MLRM: Matrix is singular!');
        print('Data Not Used!');
      }
      else if((is.positive.definite(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE ))
      {
        # print('Matrix is Positive Definite');
        lambda = solve(T_matrix, To_matrix)
        
        # reconstucting A
        if(k==1)
        {
          AUX_Pred = 
            (AUX_XLX$vectors
             %*% (lambda)
             %*% t(AUX_XLX$vectors))
        }
        else
        {
          AUX_Pred = 
            (AUX_XLX$vectors
             %*% diag(as.numeric(lambda))
             %*% t(AUX_XLX$vectors))        
        }
      }
      else
      {
        print('Matrix is not PD!');
        print('Data not used')
        AUX_Pred = Matrix::Matrix(0, N, N)
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
    
    # # Debuggin
    # print('I am here in end of Masked LRM')
    return(AUX_Pred)
}
library(igraph)
library(Matrix)
library(RSpectra)
library(dplyr)
# install.packages('lme4')
# install.packages("matrixcalc")
# library(lme4)
library(matrixcalc)
# We should find a better name for our algorithm
Masked_LRM_Algorithm <- 
  function(g, 
           CompG,
           includeSelf = 1,
           p_H = 0.1, 
           epochs = NA, 
           k = NA, 
           auxiliaryLayers = NA,
           PowerSelector = NA,
           NumberOfLeadingEigenValues = NA,
           SelectorOfEigenValues = "LA")
    {
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  E(CompG)$weight = 1
  # Network stats and adjacency matrix
  N <- vcount(g)
  AoW <- as_adjacency_matrix(g, type = "both")
  W = AoW + as_adjacency_matrix(CompG, type = "both");
  
  if(is.na(k)){
    k = N
  }
  
  if(is.na(NumberOfLeadingEigenValues)==FALSE ) {
    k = NumberOfLeadingEigenValues
  }
  
  # Allocate space for the perturbed matrix
  ScoresMatrix <- Matrix(0, N, N)
  if(includeSelf == 1)
  {
    # Call with only self to add 
    # Find the number of Eigenvalues 
    if(is.na(PowerSelector) == FALSE )
    {
      EIG_AoW <- eigen(AoW, symmetric = T, only.values = TRUE)
      EIG_AoW_Sorted = sort(abs(EIG_AoW$values), decreasing = TRUE)
      TotalEnergy = sum((EIG_AoW$values) * (EIG_AoW$values) );
      RequiredEnergy = TotalEnergy*PowerSelector;
      SumOfEnergySoFar = 0;
      k = 0
      while(SumOfEnergySoFar < RequiredEnergy)
      {
        k = k+1
        SumOfEnergySoFar = SumOfEnergySoFar + EIG_AoW_Sorted[k]*EIG_AoW_Sorted[k]
      }
    }
    
    # Diagonalise A_R
    if(is.na(k) || k == N )
    {
      XLX <- eigen(AoW, symmetric = T)
    }
    else
    {
      XLX <- eigs_sym(AoW, k = k, which = SelectorOfEigenValues)
    }
    
    # calculating lambda coefficients
    T_matrix = Matrix(0, k, k);
    To_matrix = Matrix(0, k,1);
    lambda = numeric(k);
    for(i in 1:k) {
      q_i <- XLX$vectors[, i]
      Q_I = q_i %*% t(q_i)
      To_matrix[i,1]=sum( diag( (Q_I*W)%*%(AoW) ) )
      for(j in 1:k) {
        q_j = XLX$vectors[, j]
        Q_J = q_j %*% t(q_j)
        T_matrix[i, j] = sum(diag( ( (Q_I)*W )%*%( (Q_J)*W )  ))
      }
    }
    
    if(is.singular.matrix(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE )
    {
      # TXLX = eigen(T_matrix, symmetric = T)
      print('Matrix is singular! Data not used');
      Self_score = Matrix(0, N, N)
    }
    else if(is.positive.definite(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE )
    {
      print('Matrix is Positive Definite!');
      lambda = solve(T_matrix, To_matrix)
      
      # reconstucting A
      if(k==1) 
      {
        Self_score = 
          (XLX$vectors
           %*% (lambda)
           %*% t(XLX$vectors))
      }
      else
      {
        Self_score = 
          (XLX$vectors
           %*% diag(as.numeric(lambda))
           %*% t(XLX$vectors))        
      }
    }
    else
    {
      print('Matrix is not PD! Data not used');
      Self_score = Matrix(0, N, N)
    }
    ScoresMatrix = ScoresMatrix + Self_score
  }
  
  
  
  # for each extra layer, we calculate the score
  # If there are no extra layers, the method is just SPM
  if(is.na(auxiliaryLayers) == FALSE) {
    for(iAuxL in 1:length(auxiliaryLayers)){
      AUX_Matrix <- as_adjacency_matrix(auxiliaryLayers[[iAuxL]], type = "both")
      AUX_score <- Matrix(0, N, N)
      
      # Find the number of Eigenvalues 
      if(is.na(PowerSelector) == FALSE )
      {
        EIG_A_R <- eigen(A_R, symmetric = T, only.values = TRUE)
        EIG_A_R_Sorted = sort(abs(EIG_A_R$values), decreasing = TRUE)
        TotalEnergy = sum((EIG_A_R$values) * (EIG_A_R$values) );
        RequiredEnergy = TotalEnergy*PowerSelector;
        SumOfEnergySoFar = 0;
        k = 0
        while(SumOfEnergySoFar < RequiredEnergy)
        {
          k = k+1
          SumOfEnergySoFar = SumOfEnergySoFar + EIG_A_R_Sorted[k]*EIG_A_R_Sorted[k]
        }
      }
      
      # Diagonalise AUX_score
      if(is.na(k) || k == N){
        AUX_XLX <- eigen(AUX_Matrix, symmetric = T)
      }
      else{
        AUX_XLX <- eigs_sym(AUX_Matrix, k = k, which = SelectorOfEigenValues)
      }
      
      # calculating lambda coefficients
      T_matrix = Matrix(0, k, k);
      To_matrix = Matrix(0, k,1);
      lambda = Matrix(0, k, 1);
      for(i in 1:k) {
        q_i <- AUX_XLX$vectors[, i]
        Q_I = q_i %*% t(q_i)
        To_matrix[i,1]=sum( diag( (Q_I*W)%*%(AoW) ) )
        for(j in 1:k) {
          q_j = AUX_XLX$vectors[, j]
          Q_J = q_j %*% t(q_j)
          T_matrix[i, j] = sum(diag( ( (Q_I)*W )%*%( (Q_J)*W )  ))
        }
      }
      if(is.singular.matrix(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE )
      {
        AUX_score = Matrix(0, N, N)
        print('Matrix is singular in auxillary part!');
      }
      else if((is.positive.definite(matrix(T_matrix, nrow(T_matrix), ncol(T_matrix) ) ) ==TRUE ))
      {
        print('Matrix is Positive Definite');
        lambda = solve(T_matrix, To_matrix)
        
        # reconstucting A
        if(k==1)
        {
          AUX_score = 
            (AUX_XLX$vectors
             %*% (lambda)
             %*% t(AUX_XLX$vectors))
        }
        else
        {
          AUX_score = 
            (AUX_XLX$vectors
             %*% diag(as.numeric(lambda))
             %*% t(AUX_XLX$vectors))        
        }
      }
      else
      {
        print('Matrix is not PD! Data not used');
        AUX_score = Matrix(0, N, N)
      }
      
      

      ScoresMatrix = ScoresMatrix + AUX_score
    }
  }
  
  ObservedGraph = graph_from_adjacency_matrix(W, mode = "undirected", diag = FALSE)
  non_edges <- get_non_edges(ObservedGraph)
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = ScoresMatrix[non_edges]) %>% arrange(desc(scr))
  # print(paste0('Number of Eigenvalues use is: ', toString(k))) # //??????
  return(prediction)
}
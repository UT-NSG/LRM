library(igraph)
library(Matrix)
library(RSpectra)
library(dplyr)
# We should find a better name for our algorithm
LRM_Algorithm <- 
  function(g, p_H = 0.1, epochs = 10, k = NA, 
           auxiliaryLayers = NA,
           PowerSelector = NA,
           NumberOfLeadingEigenValues = NA,
           SelectorOfEigenValues = "LA")
    {
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  L <- ecount(g)
  N <- vcount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
  if(is.na(k)){
    k = N
  }
  
  if(is.na(NumberOfLeadingEigenValues)==FALSE ) {
    k = NumberOfLeadingEigenValues
  }
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  
  # Allocate space for the perturbed matrix
  A_perturbed <- Matrix(0, N, N)
  
  for(i in 1:epochs){
    # print(paste("lp_spm", i))
    # Determine the link perturbation set dlt_E at random
    to_remove <- sample(L, links_to_remove)
    g_perturbed <- delete_edges(g, to_remove)
    
    # Compute the adjacency matrix of the graph after removal of dlt_E
    A_R <- as_adjacency_matrix(g_perturbed, type = "both")
    
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
    
    # Diagonalise A_R
    if(is.na(k) || k == N )
    {
      XLX <- eigen(A_R, symmetric = T)
    }
    else
    {
      XLX <- eigs_sym(A_R, k = k, which = SelectorOfEigenValues)
    }
    
    # Compute the adjacency matrix of the perturbation set
    dlt_A <- A - A_R
    
    # Compute correction terms for eigenvalues based on dlt_A
    dlt_l <- rowSums((t(XLX$vectors) %*% dlt_A) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # Perturb A_R but keep eigenvectors unchanged
    if(k == 1)
    {
      A_perturbed <- A_perturbed + 
        (XLX$vectors %*% (XLX$values + dlt_l) %*% t(XLX$vectors))
    }
    else
    {
      A_perturbed <- A_perturbed + 
        (XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors))
    }

  }
  # Compute the average of the perturbed matrix
  # If A is highly regular, the random removal dlt_E should not change its 
  # structure too much and A and A_perturbed should be close to each other
  A_perturbed <- A_perturbed / epochs
  # A_perturbed is the prediction score using Pure SPM method
  # Now we calculate the score 
  # using help of other auxiliary layers
  # And we simply add //????? the scores 
  # to calculate the final score of a certain link
  
  ScoresMatrix = A_perturbed;
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
      # calculating alpha_i coefficients 
      alpha <- c()
      for(i in 1:k) {
        x_i <- AUX_XLX$vectors[, i]
        alpha <- c(alpha, sum(diag((x_i %*% t(x_i)) %*% A)))
      }
      if(k==1)
      {
        AUX_score = 
          (AUX_XLX$vectors
           %*% (alpha)
           %*% t(AUX_XLX$vectors))
      }
      else
      {
        AUX_score = 
          (AUX_XLX$vectors
           %*% diag(alpha)
           %*% t(AUX_XLX$vectors))        
      }

      ScoresMatrix = ScoresMatrix + AUX_score
    }
  }
  non_edges <- get_non_edges(g)
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = ScoresMatrix[non_edges]) %>% arrange(desc(scr))
  # print(paste0('Number of Eigenvalues use is: ', toString(k))) # //??????
  return(prediction)
}
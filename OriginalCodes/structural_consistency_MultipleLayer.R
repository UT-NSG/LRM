# Adding needed libraries 
library(igraph)
library(Matrix)
library(RSpectra)
library(dplyr)

structural_consistency_MultipleLayer <- 
  function(g, p_H = 0.1, epochs = 100, k = NA, auxiliaryLayers = NA){
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  N <- vcount(g)
  L <- ecount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
  if(is.na(k)==TRUE)
  {
    k = N
  }
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  
  # Allocate space for the result
  sigma_c <- numeric(epochs)
  
  ScoresMatrix = matrix(0, N, N);
  # for each extra layer, we calculate the score
  # If there are no extra layers, the method is just SPM
  if(is.na(auxiliaryLayers) == FALSE) {
    for(iAuxL in 1:length(auxiliaryLayers)){
      TempauxiliaryLayers = auxiliaryLayers[[iAuxL]];
      AUX_Matrix <- as_adjacency_matrix(TempauxiliaryLayers[[1]], type = "both")
      AUX_score <- Matrix(0, N, N)
      # Diagonalise AUX_score
      if(is.na(k) || k == N){
        AUX_XLX <- eigen(AUX_Matrix, symmetric = T)
      }
      else{
        AUX_XLX <- eigs_sym(AUX_Matrix, k = k, which = "LA")
      }
      # calculating alpha_i coefficients 
      alpha <- c()
      for(i in 1:k) {
        x_i <- AUX_XLX$vectors[, i]
        alpha <- c(alpha, sum(diag(x_i %*% t(x_i) %*% A)))
      }
      AUX_score = 
        (AUX_XLX$vectors
         %*% diag(alpha)
         %*% t(AUX_XLX$vectors))
      ScoresMatrix = ScoresMatrix + AUX_score
    }
  }
  
  for(i in 1:epochs){
    # Determine the link perturbation set dlt_E at random
    to_remove <- sample(L, links_to_remove)
    g_perturbed <- delete_edges(g, to_remove)
    
    # Compute the adjacency matrix of the graph after removal of dlt_E
    A_R <- as_adjacency_matrix(g_perturbed, type = "both")
    
    # Diagonalise A_R
    if(is.na(k) || k == N){
      XLX <- eigen(A_R, symmetric = T)
    }else{
      XLX <- eigs_sym(A_R, k = k, which = "LA")
    }
    
    # Compute the adjacency matrix of the perturbation set
    dlt_A <- A - A_R
    
    # Compute correction terms for eigenvalues based on dlt_A
    dlt_l <- rowSums((t(XLX$vectors) %*% dlt_A) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # Perturb A_R but keep eigenvectors unchanged
    A_perturbed <- XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors)
    
    # A_perturbed = A_perturbed + ScoresMatrix
    A_perturbed = ScoresMatrix
    non_edges <- get_non_edges(g_perturbed)
    
    prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2], 
                         spm = A_perturbed[non_edges]) %>% arrange(desc(spm))
    
    # Check if the top edges are indeed part of the original graph
    prediction <- prediction[1:links_to_remove, ]
    sigma_c[i] <- sum(g[from = prediction$nodeA, to = prediction$nodeB]) /
      links_to_remove
  }
  
  return(mean(sigma_c))
}
library(igraph)
library(Matrix)
library(RSpectra)
library(dplyr)

lp_spm_with_nolinks <- function(g, p_H = 0.1, epochs = 10, k = NA, g_nolinks_known){
  k <- global_k
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  L <- ecount(g)
  N <- vcount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
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
    A_perturbed <- A_perturbed + 
      (XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors))
  }
  # Compute the average of the perturbed matrix
  # If A is highly regular, the random removal dlt_E should not change its 
  # structure too much and A and A_perturbed should be close to each other
  A_perturbed <- A_perturbed / epochs
  
  # Keep no links which are not links in both train and test data
  non_edges <- get_non_edges(g)
  g_non_edges <- graph_from_edgelist(non_edges, directed = FALSE)
  g_test_nolinks <- g_non_edges - g_nolinks_known
  test_nolinks <- as.matrix(get.data.frame(g_test_nolinks))
  # browser()
  prediction <- tibble(nodeA = test_nolinks[, 1], nodeB = test_nolinks[, 2],
                       scr = A_perturbed[test_nolinks]) %>% arrange(desc(scr))
  return(prediction)
}
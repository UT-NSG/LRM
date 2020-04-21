# Change the file name if needed
# This function gets Matrix AoM, M, and B, and returns 

library(Matrix)
library(dplyr)
library(RSpectra)

lp_reconstruct_A_using_projection_of_B_Algorithm_1_2 <- function(gT, p_H = 0.1, epochs = 10, k = NA){
  k <- global_k
  p_H_A <- global_mux_pH
  # Make weights 1 to ensure functionality
  E(gT)$weight <- 1
  # E(global_gA)$weight <- 1
  
  # Network stats and adjacency matrix
  L <- ecount(gT)
  # L_A <- ecount(global_gA)
  N <- vcount(gT)
  A <- as_adjacency_matrix(gT, type = "both")
  
  AA <- as_adjacency_matrix(global_gA, type = "both")
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  # links_to_remove_A <- round(L_A*p_H_A)
  
  # Allocate space for the perturbed matrix
  A_perturbed <- Matrix(0, N, N)
  
  for(i in 1:epochs){
    # print(paste("lp_mux_spm", i))
    # Determine the link perturbation set dlt_E at random
    to_remove <- sample(L, links_to_remove)
    g_perturbed <- delete_edges(gT, to_remove)
    # to_remove_A <- sample(L_A, links_to_remove_A)
    # gA_perturbed <- delete_edges(global_gA, to_remove_A)
    
    # Compute the adjacency matrix of the graph after removal of dlt_E
    A_R <- as_adjacency_matrix(g_perturbed, type = "both")
    
    # AA_R <- as_adjacency_matrix(gA_perturbed, type = "both")
    
    # Diagonalise A_R
    if(is.na(k) || k == N){
      XLX <- eigen(A_R, symmetric = T)
    }else{
      XLX <- eigs_sym(A_R, k = k, which = "LA")
    }
    
    # Diagonalise AA
    if(is.na(k) || k == N){
      A_XLX <- eigen(AA, symmetric = T)
    }else{
      A_XLX <- eigs_sym(AA, k = k, which = "LA")
    }
    
    # Compute the adjacency matrix of the perturbation set
    dlt_A <- A - A_R
    # dlt_AA <- AA - AA_R
    
    # Compute correction terms for eigenvalues based on dlt_A
    dlt_l <- rowSums((t(XLX$vectors) %*% dlt_A) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # coeff_calc <- function(row) {
    #   return(sum(diag(row %*% row %*% AA)))
    # }
    # print(apply(t(A_XLX$vectors), FUN = function(x) coeff_calc(x)))
    
    #############################################
    ########### Changed Parts, beside commenting lines of SPM method ####################
    AoM = A_R
    
    # creates a graph, corresponding to M
    # Since the dataset 
    #  do not distinguish between unknown and 0 edges,
    #  we have no option but to only use known edges in AoM
    #  as known edges, and create M accordingly
    ############### To change Later: when the Masking matrix is supported
    
    g_M <- g_perturbed
    M <- as_adjacency_matrix(g_M, type = "both")
    # lambda <- c()
    
    T_matrix_of_T_kj =  Matrix(0, k, k)
    To_vector_of_T_i =  Matrix(0, k, 1)
    for(i_k in 1:k) {
      q_k <- A_XLX$vectors[, i_k]
      To_vector_of_T_i[i_k] = sum(diag( ( ((q_k %*% t(q_k))*M) %*% AoM) ) ) 
      for (i_j in 1:k) {
        q_j <-A_XLX$vectors[, i_j]
        T_matrix_of_T_kj[i_k,i_j] = 
          sum(
            diag(
              (( q_k %*% t(q_k) )*M )
              %*%( ( q_j %*% t(q_j) )*M )))
      }
      
      # lambda_i = trace(x_i * x_i^T * A)
      # lambda <- c(lambda, sum(diag(x_i %*% t(x_i) %*% A)))
    }
    # if(is.positive.definite(T_matrix_of_T_kj, tol=1e-8)==FALSE){
    #   stop('T is not positive definite')
    # }
    ############ To change Later: All eigenvalues should be calculated and then the biggest k values should be chosen.
    lambda = solve(T_matrix_of_T_kj, To_vector_of_T_i)
    lambda <- as.matrix(lambda)
    ############################################
    # print(lambda)
    # print(A_XLX$values)
    
    # Compute correction terms for eigenvalues of AUX LAYER based on dlt_A
    # A_dlt_l <- rowSums((t(XLX$vectors) %*% dlt_AA) * 
    #                      t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # Perturb A_R but keep target (original) layer eigenvectors unchanged
    # A_perturbed <- A_perturbed +
    #   (XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors))
    
    # Perturb AA but keep aux layer eigenvectors unchanged
    A_perturbed <- A_perturbed + 
      (A_XLX$vectors %*% diag(c(t(lambda))) %*% t(A_XLX$vectors))
  }
  # Compute the average of the perturbed matrix
  # If A is highly regular, the random removal dlt_E should not change its 
  # structure too much and A and A_perturbed should be close to each other
  A_perturbed <- A_perturbed / epochs
  
  non_edges <- get_non_edges(gT)
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = A_perturbed[non_edges]) %>% arrange(desc(scr))
  return(prediction)
}
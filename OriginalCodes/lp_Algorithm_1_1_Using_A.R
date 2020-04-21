# Change the file name if needed
# This function gets Matrix AoM, M, and B, and returns 

library(Matrix)
library(dplyr)
library(RSpectra)

lp_Algorithm_1_1_Using_A <- function(gT, p_H = 0.1, epochs = 1, k = NA){
  # epochs = 1 because there is no sampling, and hence, there is no need to iterate upon different epochs.
  k <- global_k
  p_Mute <- global_mux_pH
  # Make weights 1 to ensure functionality
  E(gT)$weight <- 1
  
  # Network stats and adjacency matrix
  L <- ecount(gT)
  N <- vcount(gT)
  A <- as_adjacency_matrix(gT, type = "both")
  # Allocate space for the perturbed matrix
  A_perturbed <- Matrix(0, N, N)
  
  for(i in 1:epochs){
    # Diagonalise A
    if(is.na(k) || k == N){
      XLX <- eigen(A, symmetric = T)
    }else{
      XLX <- eigs_sym(A, k = k, which = "LA")
    }
    
    
    #############################################
    ########### Changed Parts, beside commenting lines of SPM method ####################
    AoM = A
    
    # creates a graph, corresponding to M
    # Since the dataset 
    #  do not distinguish between unknown and 0 edges,
    #  we have no option but to only use known edges in AoM
    #  as known edges, and create M accordingly
    ############### To change Later: when the Masking matrix is supported
    
    g_M <- gT
    M <- as_adjacency_matrix(g_M, type = "both")
    # lambda <- c()
    
    T_matrix_of_T_kj =  Matrix(0, k, k)
    To_vector_of_T_i =  Matrix(0, k, 1)
    for(i_k in 1:k) {
      q_k <- XLX$vectors[, i_k]
      To_vector_of_T_i[i_k] = sum(diag( ( ((q_k %*% t(q_k))*M) %*% AoM) ) ) 
      for (i_j in 1:k) {
        q_j <-XLX$vectors[, i_j]
        T_matrix_of_T_kj[i_k,i_j] = 
          sum(
            diag(
              (( q_k %*% t(q_k) )*M )
              %*%( ( q_j %*% t(q_j) )*M )))
      }
    }
    # if(is.positive.definite(T_matrix_of_T_kj, tol=1e-8)==FALSE){
    #   stop('T is not positive definite')
    # }
    ############ To change Later: All eigenvalues should be calculated and then the biggest k values should be chosen.
    lambda = solve(T_matrix_of_T_kj, To_vector_of_T_i)
    lambda <- as.matrix(lambda)
    ############################################
    # print(lambda)
    
    # Perturb A but keep aux layer eigenvectors unchanged
    A_perturbed <- A_perturbed + 
      (XLX$vectors %*% diag(c(t(lambda))) %*% t(XLX$vectors))
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
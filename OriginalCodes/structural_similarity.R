# Measurment of Structural Similarity of layers

structural_similarity <- function(gT, k = NA){
 
   # Network stats and adjacency matrix
  N <- vcount(gT)
  L <- ecount(gT)
  L_A <- ecount(global_gA)
  
  A <- as_adjacency_matrix(gT, type = "both")
  B <- as_adjacency_matrix(global_gA, type = "both")
  
  #k <- 10
  
  # Diagonalise A
  if(is.na(k) || k == N){
    A_XLX <- eigen(A, symmetric = T)
  }else{
    A_XLX <- eigs_sym(A, k = k, which = "LA")
  }
  
  # Diagonalise B
  if(is.na(k) || k == N){
    B_XLX <- eigen(B, symmetric = T)
  }else{
    B_XLX <- eigs_sym(B, k = k, which = "LA")
  }
  
  SS <- abs(t(A_XLX$vectors) %*% B_XLX$vectors)
  
  return(SS)
}

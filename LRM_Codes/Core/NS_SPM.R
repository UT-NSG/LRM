# for E, ecount, vcount, as_adjacency_matrix, laplacian_matrix, and delete_edges
library(igraph)

# for Matrix
library(Matrix)

# for eigs_sym
library(RSpectra)


# library(dplyr)

NS_SPM = 
  function(
    gLT,
    gNLT,
    auxL,
    NumberOfLeadingEigenValues,
    SelectorOfEigenValues,
    PowerSelector,
    EPType = 'None',
    # Local variables to SPM
    p_H = 0.1,
    local_epochs = 10
    )
    {
    
    # rint(auxL)
    if(auxL=='Train') {
      g = gLT
    }
    else {
      g = auxL
    }
    
    # Make weights 1 to ensure functionality
    igraph::E(g)$weight = 1
    
    # Network stats and adjacency matrix
    L = igraph::ecount(g)
    N = igraph::vcount(g)
    
    
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
    
    AUX = matrix_creator_Function_Hanlge(g)
  
    if(is.na(NumberOfLeadingEigenValues)==FALSE ) {
      k = NumberOfLeadingEigenValues
    }
    else if(is.na(NumberOfLeadingEigenValues)==TRUE) {
      k = N
    }
    
    # Size of the probe set E^p
    links_to_remove = round(L*p_H)
    
    # Allocate space for the perturbed matrix
    AUX_perturbed = Matrix::Matrix(0, N, N)
    
    for(i in 1:local_epochs){
      # Determine the link perturbation set dlt_E at random
      to_remove = sample(L, links_to_remove)
      gAUX_perturbed = igraph::delete_edges(g, to_remove)
      
      # make sure the size stays the same 
      if(igraph::vcount(gAUX_perturbed) < N) {
        gAUX_perturbed = 
          igraph::add_vertices(gAUX_perturbed, N - igraph::vcount(gAUX_perturbed))
      }
      
      # Compute the adjacency matrix of the graph after removal of dlt_E
      AUX_R = matrix_creator_Function_Hanlge(gAUX_perturbed)
      
      # Find the number of Eigenvalues 
      if(is.na(PowerSelector) == FALSE )
      {
        EIG_AUX_R = eigen(AUX_R, symmetric = T, only.values = TRUE)
        EIG_AUX_R_Sorted = sort(abs(EIG_AUX_R$values), decreasing = TRUE)
        TotalEnergy = sum((EIG_AUX_R$values) * (EIG_AUX_R$values) );
        RequiredEnergy = TotalEnergy*PowerSelector;
        SumOfEnergySoFar = 0;
        k = 0
        while(SumOfEnergySoFar < RequiredEnergy)
        {
          k = k+1
          SumOfEnergySoFar = SumOfEnergySoFar + EIG_AUX_R_Sorted[k]*EIG_AUX_R_Sorted[k]
        }
      }
      
      # Diagonalise AUX_R
      if(k == N )
      {
        XLX = eigen(AUX_R, symmetric = T)
      }
      else
      {
        XLX = RSpectra::eigs_sym(AUX_R, k = k, which = SelectorOfEigenValues)
      }
      
      # Compute the adjacency matrix of the perturbation set
      dlt_AUX = AUX - AUX_R
      
      # Compute correction terms for eigenvalues based on dlt_A
      dlt_l = 
        rowSums((t(XLX$vectors) %*% dlt_AUX) * 
                         t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
      
      # print(XLX$values + dlt_l)
      # print(dlt_l)
      
      # Perturb AUX_R but keep eigenvectors unchanged
      if(k == 1)
      {
        AUX_perturbed = 
          AUX_perturbed + 
          (XLX$vectors %*% (XLX$values + dlt_l) %*% t(XLX$vectors))
      }
      else
      {
        AUX_perturbed = 
          AUX_perturbed + 
          (XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors))
      }
    }
    # Compute the average of the perturbed matrix
    # If AUX is highly regular, the random removal dlt_E should not change its 
    # structure too much and A and A_perturbed should be close to each other
    AUX_predection = AUX_perturbed / local_epochs
    
    # print(paste0('Number of Eigenvalues used is: ', toString(k))) # //??????
    return(AUX_predection)
}
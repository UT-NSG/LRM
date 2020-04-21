library(purrr)
library(LinkPrediction)
source("prune_predict_assess_MultiLayer.R")

prune_recover_MultipleLayer <- 
  function (g, ..., 
            probes = seq(0.1, 0.9, 0.1), 
            epochs = 30, 
            preserve_conn = FALSE, 
            use_weights = FALSE, 
            auxiliaryLayers = NA,
            PowerSelector = NA,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA",
            includeOnlyLRM = 0) 
{
  predictors <- list(...)
  if (("lp_mce" %in% predictors || "lp_leig" %in% predictors || 
       "lp_isomap" %in% predictors) && !preserve_conn) {
    stop(paste0("One or more of your link predictors requires a connected ", 
                "network, set preserve_conn to TRUE!"))
  }
  if (preserve_conn) {
    if (use_weights) {
      g_mst <- mst(g)
    }
    else {
      g_mst <- mst(g, algorithm = "unweighted")
    }
    removable_links <- E(g)[-g[from = head_of(g_mst, E(g_mst)), 
                               to = tail_of(g_mst, E(g_mst)), edges = TRUE]]
  }
  else {
    removable_links <- E(g)
  }
  col_len <- length(predictors) * length(probes) * epochs
  # Creating the result of assesment tibble, to store the test result in it
  res <- 
    tibble(
      method = rep(unlist(predictors), length(probes) * epochs),
      epoch = rep(1:epochs, each = length(predictors) * length(probes)),
      frac_rem = rep(rep(probes, each = length(predictors)), epochs), 
      links_rem = round(frac_rem * length(removable_links)), 
      recall_at_k = numeric(col_len), 
      aupr = numeric(col_len), 
      auroc = numeric(col_len), 
      avg_prec = numeric(col_len)
      )
  
  evaluation <- map2(as.list(res$method), 
                     as.list(res$links_rem), 
                     prune_predict_assess_MultiLayer, 
                     g = g, 
                     removable_links = removable_links,
                     auxiliaryLayers=auxiliaryLayers,
                     PowerSelector = PowerSelector,
                     NumberOfLeadingEigenValues = NumberOfLeadingEigenValues,
                     SelectorOfEigenValues = SelectorOfEigenValues,
                     includeOnlyLRM = includeOnlyLRM
                     )
  
  res$recall_at_k <- map_dbl(evaluation, function(x) x["recall_at_k"])
  res$aupr <- map_dbl(evaluation, function(x) x["aupr"])
  res$auroc <- map_dbl(evaluation, function(x) x["auroc"])
  res$avg_prec <- map_dbl(evaluation, function(x) x["avg_prec"])
  return(res)
}

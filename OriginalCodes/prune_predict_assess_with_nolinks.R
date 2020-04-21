library(precrec)
prune_predict_assess_with_nolinks <- function (g, method, links_rem, removable_links) 
{
  perf <- numeric(4)
  names(perf) <- c("recall_at_k", "aupr", "auroc", "avg_prec")
  edg <- sample(removable_links, links_rem)
  g_perturbed <- delete_edges(g, edg)
  
  g_nolinks <- graph_from_edgelist(get_non_edges(g), directed = FALSE)
  edg_nolinks <- sample(E(g_nolinks), links_rem)
  N_perturbed <- delete_edges(g_nolinks, edg_nolinks)
  
  pred <- do.call(method, list(g = g_perturbed, g_nolinks_known = N_perturbed))
  scr <- nrow(pred):1
  classes <- g[from = pred$nodeA, to = pred$nodeB]
  perf["recall_at_k"] <- sum(classes[1:links_rem])/links_rem
  meval <- evalmod(scores = scr, labels = classes)
  perf["aupr"] <- attr(meval$prcs[[1]], "auc")
  perf["auroc"] <- attr(meval$rocs[[1]], "auc")
  perf["avg_prec"] <- mean(meval$prcs[[1]]$y[meval$prcs[[1]]$orig_points])
  return(perf)
}

library(precrec)
prune_predict_assess_MultiLayer_Masked_LRM_Compatible <- 
  function (g, 
            method, 
            links_rem, 
            removable_links, 
            NoLinksToLinksRatio = 2,
            includeSelf = 1,
            auxiliaryLayers = NA, 
            PowerSelector = NA,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA") 
{
  perf <- numeric(4)
  names(perf) <- c("recall_at_k", "aupr", "auroc", "avg_prec")
  edg <- sample(removable_links, links_rem)
  g_train_Links <- delete_edges(g, edg)
  
  print(paste0("Number of Removed Links: ", links_rem) )
  
  CompG = complementer(g, loops = FALSE)
  removable_Nolinks = E(CompG)
  Noedg <- sample(removable_Nolinks, links_rem*NoLinksToLinksRatio)
  g_train_NoLinks <- delete_edges(CompG, Noedg)
  
  pred <- 
    do.call(
      method, 
      list(
        g = g_train_Links,
        CompG = g_train_NoLinks,
        includeSelf = includeSelf,
        auxiliaryLayers = auxiliaryLayers, 
        NumberOfLeadingEigenValues = NumberOfLeadingEigenValues, 
        PowerSelector = PowerSelector,
        SelectorOfEigenValues = SelectorOfEigenValues
        )
      )
  scr <- nrow(pred):1
  classes <- g[from = pred$nodeA, to = pred$nodeB]
  perf["recall_at_k"] <- sum(classes[1:links_rem])/links_rem
  meval <- evalmod(scores = scr, labels = classes)
  perf["aupr"] <- attr(meval$prcs[[1]], "auc")
  perf["auroc"] <- attr(meval$rocs[[1]], "auc")
  perf["avg_prec"] <- mean(meval$prcs[[1]]$y[meval$prcs[[1]]$orig_points])
  return(perf)
}

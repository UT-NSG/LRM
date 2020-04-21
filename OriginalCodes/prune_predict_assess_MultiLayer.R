library(precrec)
prune_predict_assess_MultiLayer <- 
  function (g, 
            method, 
            links_rem, 
            removable_links, 
            auxiliaryLayers = NA, 
            PowerSelector = NA,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA",
            includeOnlyLRM = 0) 
{
  perf <- numeric(4)
  names(perf) <- c("recall_at_k", "aupr", "auroc", "avg_prec")
  edg <- sample(removable_links, links_rem)
  g_train <- delete_edges(g, edg)
  if(method == "LinkPrediction_Algorithm0_OnlyLRMSupport"){
  pred <- 
    do.call(
      method, 
      list(
        g = g_train, 
        auxiliaryLayers = auxiliaryLayers, 
        NumberOfLeadingEigenValues = NumberOfLeadingEigenValues, 
        PowerSelector = PowerSelector,
        SelectorOfEigenValues = SelectorOfEigenValues, 
        includeOnlyLRM = includeOnlyLRM
        )
      )
  }
  else{
    pred <- 
      do.call(
        method, 
        list(
          g = g_train, 
          auxiliaryLayers = auxiliaryLayers, 
          NumberOfLeadingEigenValues = NumberOfLeadingEigenValues, 
          PowerSelector = PowerSelector,
          SelectorOfEigenValues = SelectorOfEigenValues
        )
      )
  }
  scr <- nrow(pred):1
  classes <- g[from = pred$nodeA, to = pred$nodeB]
  perf["recall_at_k"] <- sum(classes[1:links_rem])/links_rem
  meval <- evalmod(scores = scr, labels = classes)
  perf["aupr"] <- attr(meval$prcs[[1]], "auc")
  perf["auroc"] <- attr(meval$rocs[[1]], "auc")
  perf["avg_prec"] <- mean(meval$prcs[[1]]$y[meval$prcs[[1]]$orig_points])
  return(perf)
}

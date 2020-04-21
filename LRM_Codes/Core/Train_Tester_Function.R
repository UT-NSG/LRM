# for evalmod function
library(precrec)

# for delete_edges, complementer, E, and vcount
library(igraph)

# for calling the method over the parameters and aggregate the results
source('Core/Extended_Method_Caller.R')

Train_Tester_Function = 
  function (links_rem,
            target_Graph, 
            methodList,
            sourceList,
            aggregation_method = 'SimpleAddition',
            extra_processing = 'None',
            removable_links, 
            NoLinksToLinksRatio = 2,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA",
            PowerSelector = NA,
            AddingAUXAdjacencyScale = 0
            ) 
{
  perf = numeric(4)
  names(perf) = c("recall_at_k", "aupr", "auroc", "avg_prec")
  
  # print('Beginning of Train_Tester')
  
  N = igraph::vcount(target_Graph)
  
  # Creating train link set
  edg = sample(removable_links, links_rem)
  g_train_Links = igraph::delete_edges(target_Graph, edg)
  
  # make sure the size stays the same 
  if(igraph::vcount(g_train_Links) < N) {
    g_train_Links = 
      igraph::add_vertices(g_train_Links, N - igraph::vcount(g_train_Links))
  }
  
  
  # Creating train nolink set
  CompG = igraph::complementer(target_Graph, loops = FALSE)
  removable_Nolinks = igraph::E(CompG)
  NumberOfNoLinksToRemove = min(links_rem*NoLinksToLinksRatio, igraph::vcount(CompG))
  # print('Number of no links to remove')
  # print(NumberOfNoLinksToRemove)
  # print('Number of edges in G complement')
  # print(igraph::vcount(CompG))
  Noedg = 
    sample(
      removable_Nolinks, 
      NumberOfNoLinksToRemove
    )
  # if(NumberOfNoLinksToRemove==igraph::vcount(CompG)) {
    # print('All non-edges were used.')
  # }
  g_train_NoLinks = igraph::delete_edges(CompG, Noedg)
  
  # make sure the size stays the same 
  if(igraph::vcount(g_train_NoLinks) < N) {
    g_train_NoLinks = 
      igraph::add_vertices(g_train_NoLinks, N - igraph::vcount(g_train_NoLinks))
  }
  
  # print('Before Extended Method Caller of Train_Tester')
  # Calling the actual score aggregator
  pred = Extended_Method_Caller(
    g_train_Links = g_train_Links, 
    g_train_NoLinks = g_train_NoLinks,
    methodList = methodList,
    sourceList = sourceList,
    aggregation_method = aggregation_method,
    extra_processing = extra_processing,
    NumberOfLeadingEigenValues = NumberOfLeadingEigenValues, 
    SelectorOfEigenValues = SelectorOfEigenValues,
    PowerSelector = PowerSelector,
    AddingAUXAdjacencyScale = AddingAUXAdjacencyScale
  )
  
  # # Debuggin
  # print('I am here after call to Extended_Method_Caller in Train_Tester')
  
  # ?????? # what is the next line doing? Is it only using the order, instead of the actual score calculated by the method?
  scr = nrow(pred):1
  classes = target_Graph[from = pred$nodeA, to = pred$nodeB]
  perf["recall_at_k"] = sum(classes[1:links_rem])/links_rem
  meval = precrec::evalmod(scores = scr, labels = classes)
  perf["aupr"] = attr(meval$prcs[[1]], "auc")
  perf["auroc"] = attr(meval$rocs[[1]], "auc")
  perf["avg_prec"] = mean(meval$prcs[[1]]$y[meval$prcs[[1]]$orig_points])
  return(perf)
}

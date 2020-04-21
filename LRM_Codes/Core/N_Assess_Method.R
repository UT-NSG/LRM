# In the Name of God

#### Preparations ####
## Adding needed libraries 
# for vcount and E function
library(igraph)

# for map function
library(purrr)

# library(LinkPrediction)
# library(lattice)

## Adding needed source codes 
# Creates train and test sets, pass the train to the method caller and calculate accuracy based on test set.
source('Core/Train_Tester_Function.R')

# It is advised to pass same size graphes to the function
N_Assess_Method = 
  function (Tg_Graph, 
            methodList,
            sourceList, 
            aggregation_method = 'SimpleAddition',
            extra_processing = 'None',
            legendFig,
            probesOfLinks = seq(0.1, 0.9, 0.1), 
            epochs = 30,
            NumberOfLeadingEigenValues = NA,
            SelectorOfEigenValues = "LA",
            PowerSelector = NA,
            NoLinksToLinksRatio=2,
            AddingAUXAdjacencyScale = 0)
    {
    
    # print('I am at the beginning of Assess_Method')
    # Making sure the graphs have the same size
    NumSourceLay = length(sourceList)
    
    # Making sure that all graphes are the same size 
    maximumVcount = igraph::vcount(Tg_Graph)
    for (i_m_c in 1:NumSourceLay) {
      if(sourceList[i_m_c]=='Train') {
        # noting is needed to be done
      }
      else {
        maximumVcount = 
          max(maximumVcount, igraph::vcount(sourceList[[i_m_c]]) )
      }
    }
    if(igraph::vcount(Tg_Graph) < maximumVcount) {
      Tg_Graph = 
        igraph::add_vertices(Tg_Graph, maximumVcount - igraph::vcount(Tg_Graph))
    }
    for (i_m_v in 1:NumSourceLay) {
      if(sourceList[[i_m_v]]=='Train') {
        # noting to do
      }
      else {
        if(igraph::vcount(sourceList[[i_m_v]]) < maximumVcount) {
          sourceList[[i_m_v]] = 
            igraph::add_vertices(sourceList[[i_m_v]], 
                         maximumVcount - igraph::vcount(sourceList[[i_m_v]]))
        }
      }
    }
    
    #Creating the removable links list
    removable_links = igraph::E(Tg_Graph)
    # finding the length of result tibble
    col_len = length(probesOfLinks) * epochs
    
    # print('I am before assesment in Assess_Method')
    # Creating the result of assesment tibble, to store the test result in it
    assessment = 
      tibble::tibble(
        method = rep(legendFig, length(probesOfLinks) * epochs),
        epoch = rep(1:epochs, each = length(probesOfLinks)),
        frac_rem = rep(probesOfLinks, epochs), 
        links_rem = round(frac_rem * length(removable_links)), 
        NoLToLRatio = rep(NoLinksToLinksRatio,times=length(epoch)),
        recall_at_k = numeric(col_len), 
        aupr = numeric(col_len), 
        auroc = numeric(col_len), 
        avg_prec = numeric(col_len)
      )

    # print('I am after  assesment in Assess_Method')
    
    evaluation = map(
        as.list(assessment$links_rem), 
        Train_Tester_Function, 
        target_Graph = Tg_Graph, 
        methodList = methodList,
        sourceList = sourceList,
        aggregation_method = aggregation_method,
        extra_processing = extra_processing,
        removable_links = removable_links,
        NoLinksToLinksRatio = NoLinksToLinksRatio,
        NumberOfLeadingEigenValues = NumberOfLeadingEigenValues,
        SelectorOfEigenValues = SelectorOfEigenValues,
        PowerSelector = PowerSelector,
        AddingAUXAdjacencyScale = AddingAUXAdjacencyScale
    )
    
    # print('I am after  evaluations in Assess_Method')
    assessment$recall_at_k = 
      purrr::map_dbl(evaluation, function(x) x["recall_at_k"])
    assessment$aupr = 
      purrr::map_dbl(evaluation, function(x) x["aupr"])
    assessment$auroc =
      purrr::map_dbl(evaluation, function(x) x["auroc"])
    assessment$avg_prec =
      purrr::map_dbl(evaluation, function(x) x["avg_prec"])
    
    # print('I am at the end of Assess_Method')
    return(assessment)
}


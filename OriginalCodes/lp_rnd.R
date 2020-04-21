lp_rnd <- function(g){
  non_edges <- get_non_edges(g)
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = sample(nrow(non_edges))) %>% 
    arrange(desc(scr))
  return(prediction)
}
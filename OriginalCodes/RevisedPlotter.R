library(ggplot2)
library(dplyr) 


RevisedPlotter <- function(res, metric = "recall_at_k", colours = NA, 
                                err = "none", Version = "FracRem"){
  
  if( Version == "FracRem"){
    y_axis <- switch(metric, 
                     recall_at_k = "Recall@k", 
                     aupr = "AUPR", 
                     auroc = "AUROC", 
                     avg_prec = "Avg. Precision",
                     stop("The selected metric is not valid!")
    )
    if(any(is.na(colours))){
      colours = seq_along(unique(res$method))
    }
    else if(length(colours) != length(unique(res$method))){
      stop("The number of link prediction methods and colours does not match!!!")
    }
    if(!(err %in% c("none", "sd", "se"))){
      stop("The indicated error type is not valid!")
    }
    
    res <- res %>% mutate(perf = res[[metric]])
    
    if(err == "none"){
      res <- res %>% group_by(method, frac_rem) %>% 
        summarise(avg_perf = mean(perf))
      err_bar_width <- 0
    }
    else if(err == "sd"){
      res <- res %>% group_by(method, frac_rem) %>% 
        summarise(avg_perf = mean(perf), 
                  up_err = avg_perf + sd(perf), 
                  dw_err = avg_perf - sd(perf))
      err_bar_width <- 0.01
    }
    else if(err == "se"){
      res <- res %>% group_by(method, frac_rem) %>% 
        summarise(avg_perf = mean(perf), 
                  up_err = avg_perf + sd(perf)/sqrt(n()), 
                  dw_err = avg_perf - sd(perf)/sqrt(n()))
      err_bar_width <- 0.01
    }
    
    p <- res %>% ggplot(aes_(~frac_rem, ~avg_perf, colour = ~method)) + 
      geom_line() + 
      geom_errorbar(aes_(ymin = ~dw_err, ymax = ~up_err), width = err_bar_width) +
      scale_colour_manual(values = colours) +
      labs(x = "Fraction of removed links", 
           y = y_axis, colour = "Predictor") + 
      theme_bw() + theme(legend.background = element_blank(), 
                         legend.position = "top")
  }
  
  if(Version == "NumEig"){
    y_axis <- switch(metric, 
                     recall_at_k = "Recall@k", 
                     aupr = "AUPR", 
                     auroc = "AUROC", 
                     avg_prec = "Avg. Precision",
                     stop("The selected metric is not valid!")
    )
    if(any(is.na(colours))){
      colours = seq_along(unique(res$method))
    }
    else if(length(colours) != length(unique(res$method))){
      stop("The number of link prediction methods and colours does not match!!!")
    }
    if(!(err %in% c("none", "sd", "se"))){
      stop("The indicated error type is not valid!")
    }
    
    res <- res %>% mutate(perf = res[[metric]])
    
    if(err == "none"){
      res <- res %>% group_by(method, k) %>% 
        summarise(avg_perf = mean(perf))
      
      err_bar_width <- 0
    }
    else if(err == "sd"){
      res <- res %>% group_by(method, k) %>% 
        summarise(avg_perf = mean(perf), 
                  up_err = avg_perf + sd(perf), 
                  dw_err = avg_perf - sd(perf))
      err_bar_width <- 0.01
    }
    else if(err == "se"){
      res <- res %>% group_by(method, k) %>% 
        summarise(avg_perf = mean(perf), 
                  up_err = avg_perf + sd(perf)/sqrt(n()), 
                  dw_err = avg_perf - sd(perf)/sqrt(n()))
      err_bar_width <- 0.01
    }
    
    p <- res %>% ggplot(aes_(~k, ~avg_perf, colour = ~method)) + 
      geom_line() + 
      geom_errorbar(aes_(ymin = ~dw_err, ymax = ~up_err), width = err_bar_width) +
      scale_colour_manual(values = colours) +
      labs(x = "Number of Leading Eigenvalues", 
           y = y_axis, colour = "Predictor") + 
      theme_bw() + theme(legend.background = element_blank(), 
                         legend.position = "top")
  }
  
  
  return(p)
}
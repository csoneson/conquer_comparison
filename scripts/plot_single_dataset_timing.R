plot_timing <- function(timinglist, colvec, exts, summary_data = list()) {
  ## Extract elapsed time for each method and data set instance
  timings <- sapply(timinglist, function(i) {
    i["user.self"] + i["sys.self"] + i["user.child"] + i["sys.child"]
  })
  
  ## Create data frame with all timing information
  timings_full <- data.frame(method = names(timings), timing = timings) %>%
    tidyr::separate(method, into = c("method", "ncells", "repl", "user", "self"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(exts, "", method))
  
  ## Plot timing as function of number of cells per group
  print(timings_full %>%    
          ggplot(aes(x = ncells, y = timing, group = method, color = method)) + 
          geom_point(alpha = 0.5) + geom_smooth(se = FALSE) + scale_y_log10() + 
          theme_bw() + xlab("Number of cells per group") + 
          scale_color_manual(values = structure(colvec, names = gsub(exts, "", names(colvec)))))
  
  summary_data$timing_full <- timings_full
  return(invisible(summary_data))
}

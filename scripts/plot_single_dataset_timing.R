plot_timing <- function(timinglist, colvec, exts = exts, summary_data = list()) {
  timings <- sapply(timinglist, function(i) i["elapsed"])
  
  timings_full <- data.frame(method = names(timings), timing = timings) %>%
    tidyr::separate(method, into = c("method", "ncells", "repl", "elapsed"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(exts, "", method))
  
  timings <- data.frame(method = names(timings), timing = timings) %>% 
    tidyr::separate(method, into = c("method", "ncells", "repl", "elapsed"), sep = "\\.") %>%
    group_by(method, ncells) %>% dplyr::summarise(timing = median(timing)) %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(method = gsub(exts, "", method))
  
  summary_data$timing <- rbind(summary_data$timing, timings)
  summary_data$timing_full <- rbind(summary_data$timing_full, timings_full)
  
  print(timings %>%    
          ggplot(aes(x = ncells, y = timing, group = method, color = method)) + 
          geom_line(size = 2.5) + scale_y_log10() + theme_bw() + 
          xlab("Number of cells per group") + 
          scale_color_manual(values = structure(colvec, names = gsub(exts, "", names(colvec)))))
  return(invisible(summary_data))
}

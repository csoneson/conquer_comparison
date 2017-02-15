source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Plot time used for each method
#' 
plot_timing <- function(timinglist, colvec, summary_data = list()) {
  timings <- sapply(timinglist, function(i) i["elapsed"])
  timings <- data.frame(method = names(timings), timing = timings) %>% 
    tidyr::separate(method, into = c("method", "nsamples", "repl", "elapsed"), sep = "\\.") %>%
    group_by(method, nsamples) %>% dplyr::summarise(timing = median(timing)) %>%
    dplyr::mutate(nsamples = as.numeric(as.character(nsamples)))
  summary_data$timing <- rbind(summary_data$timing, timings)
  print(timings %>%    
          ggplot(aes(x = nsamples, y = timing, group = method, color = method)) + 
          geom_line(size = 2.5) + scale_y_log10() + theme_bw() + 
          xlab("Number of cells") + 
          scale_color_manual(values = colvec))
  return(invisible(summary_data))
}

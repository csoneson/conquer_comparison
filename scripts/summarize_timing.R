summarize_timing <- function(figdir, datasets, exts) {
  
  ## ------------------------------- Timing ----------------------------------- ##
  pdf(paste0(figdir, "/summary_timing", exts, ".pdf"), width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/timing/", ds, exts, 
                   "_timing_summary_data.rds"))
  })
  y <- lapply(summary_data_list, function(m) {
    if (!is.null(m$timing)) {
      m$timing %>% group_by(dataset, filt, nsamples) %>%
        dplyr::mutate(timing = timing/max(timing))
    } else {
      NULL
    }
  })
  y <- do.call(rbind, y)
  ## Boxplots
  print(ggplot(y, aes(x = method, y = timing, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative timing") + 
          scale_color_manual(values = cols, name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  print(ggplot(y, aes(x = method, y = timing, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative timing") + 
          scale_color_manual(values = cols, name = "") + scale_y_log10() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  ## Barplots
  y %>% group_by(method) %>% dplyr::summarize(mean = mean(timing), sd = sd(timing)) %>%
    ggplot(aes(x = method, y = mean, fill = method)) + 
    geom_errorbar(aes(ymin = min(mean)/2, ymax = mean + sd), width = 0.2) + 
    geom_bar(stat = "identity") + 
    theme_bw() + xlab("") + ylab("Relative timing") + 
    scale_fill_manual(values = cols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ## Dependence on number of samples
  y2 <- lapply(summary_data_list, function(m) {
    if (!is.null(m$timing)) {
      m$timing
    } else {
      NULL
    }
  })
  y2 <- do.call(rbind, y2)
  print(y2 %>% group_by(method, dataset, filt) %>%
          arrange(nsamples) %>% 
          dplyr::mutate(dt = c(0, diff(timing)),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, diff(nsamples)),
                        s = c(0, nsamples[1:(length(nsamples) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = cols, name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  print(y2 %>% group_by(method, dataset, filt) %>%
          arrange(nsamples) %>% 
          dplyr::mutate(dt = c(0, timing[2:(length(timing))]),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, nsamples[2:(length(nsamples))]),
                        s = c(0, nsamples[1:(length(nsamples) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = cols, name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  dev.off()
  
  
}
summarize_timing <- function(figdir, datasets, exts) {
  
  ## ------------------------------- Timing ----------------------------------- ##
  pdf(paste0(figdir, "/summary_timing", exts, ".pdf"), width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/timing/", ds, exts, 
                   "_timing_summary_data.rds"))
  })
  y <- lapply(summary_data_list, function(m) {
    if (!is.null(m$timing)) {
      m$timing %>% group_by(dataset, filt, ncells) %>%
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
          arrange(ncells) %>% 
          dplyr::mutate(dt = c(0, diff(timing)),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, diff(ncells)),
                        s = c(0, ncells[1:(length(ncells) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = cols, name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  print(y2 %>% group_by(method, dataset, filt) %>%
          arrange(ncells) %>% 
          dplyr::mutate(dt = c(0, timing[2:(length(timing))]),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, ncells[2:(length(ncells))]),
                        s = c(0, ncells[1:(length(ncells) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = 0) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = cols, name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  ## Dependence on number of genes
  ngenes <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, 
                   "_ngenes.rds"))
  })
  ngenes <- do.call(rbind, ngenes) %>% group_by(dataset, filt, ncells) %>%
    summarize(ngenes = median(ngenes)) %>% 
    mutate(ncells = as.numeric(ncells))
  y3 <- inner_join(y2, ngenes)
  
  ns <- y3 %>% filter(method == y3$method[1]) %>% group_by(ncells) %>% tally()
  ns_keep <- ns$ncells[ns$n > 1]
  
  y3 %>% filter(ncells %in% ns_keep) %>% group_by(method, ncells, filt) %>%
    arrange(ngenes) %>%
    dplyr::mutate(dt = c(0, timing[2:(length(timing))]),
                  t = c(0, timing[1:(length(timing) - 1)]),
                  dg = c(0, ngenes[2:(length(ngenes))]),
                  g = c(0, ngenes[1:(length(ngenes)) - 1])) %>%
    dplyr::mutate(reltime = (dt/t)/(dg/g)) %>%
    filter(!is.na(reltime)) %>%
    ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = 0) + 
    geom_point(position = position_jitter(width = 0.2)) + 
    theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of genes") + 
    scale_color_manual(values = cols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ## 2d heatmap
  ngenes <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, 
                   "_ngenes.rds"))
  })
  ngenes <- do.call(rbind, ngenes) %>% ungroup() %>% mutate(ncells = as.numeric(ncells)) %>%
    mutate(repl = as.numeric(repl))
  y4 <- lapply(summary_data_list, function(m) {
   m$timing_full
  })
  y4 <- do.call(rbind, y4) %>% mutate(repl = as.numeric(repl))
  y4 <- inner_join(y4, ngenes)
  for (m in unique(y4$method)) {
    scatterplot3d(subset(y4, method == m)$ncells, 
                  subset(y4, method == m)$ngenes, 
                  subset(y4, method == m)$timing, type = "h", pch = 19)
  }
  
  dev.off()
  
}
summarize_relfprtpr <- function(figdir, datasets, exts, dtpext, cols = cols) {
  ## ---------------------------------- True FPR ------------------------------ ##
  pdf(paste0(figdir, "/summary_heatmaps_relfprtpr", exts, dtpext, ".pdf"),
      width = 10, height = 4 * length(datasets))
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/results_relativetruth/", ds, exts, "_results_relativetruth_summary_data.rds"))
  })
  ## Heatmap of relative FPRs
  y <- lapply(summary_data_list, function(m) {
    m$fpr_relative %>% 
      tidyr::separate(basemethod, c("method", "n_samples", "repl"), sep = "\\.") %>%
      dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", n_samples, ".", repl)) %>%
      dplyr::select(method, dataset, FPR)
  })
  y <- do.call(rbind, y) %>% dcast(dataset ~ method, value.var = "FPR") %>%
    tidyr::separate(dataset, c("ds", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::arrange(ds, as.numeric(as.character(n_samples))) %>% 
    dplyr::select(-ds, -filt, -n_samples, -repl)  %>% as.data.frame()
  rownames(y) <- y$dataset
  y$dataset <- NULL
  
  ## Remove extension from method name
  colnames(y) <- gsub(exts, "", colnames(y))
  
  annotation_row = data.frame(id = rownames(y)) %>% 
    tidyr::separate(id, c("dataset", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(n_samples = factor(n_samples, 
                                     levels = as.character(sort(unique(as.numeric(as.character(n_samples)))))))
  rownames(annotation_row) <- annotation_row$id
  
  pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", main = "Relative FPR",
           display_numbers = TRUE, colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = seq(0, 0.6, length.out = 101), 
           annotation_row = dplyr::select(annotation_row, n_samples, dataset), show_rownames = FALSE,
           annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
           annotation_colors = list(method = structure(cols, names = gsub(exts, "", names(cols)))[colnames(y)]),
           annotation_names_col = FALSE)
  
  ## Heatmap of relative TPRs
  y <- lapply(summary_data_list, function(m) {
    m$tpr_relative %>% 
      tidyr::separate(basemethod, c("method", "n_samples", "repl"), sep = "\\.") %>%
      dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", n_samples, ".", repl)) %>%
      dplyr::select(method, dataset, TPR)
  })
  y <- do.call(rbind, y) %>% dcast(dataset ~ method, value.var = "TPR") %>%
    tidyr::separate(dataset, c("ds", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::arrange(ds, as.numeric(as.character(n_samples))) %>% 
    dplyr::select(-ds, -filt, -n_samples, -repl)  %>% as.data.frame()
  rownames(y) <- y$dataset
  y$dataset <- NULL
  
  ## Remove extension from method name
  colnames(y) <- gsub(exts, "", colnames(y))
  
  annotation_row = data.frame(id = rownames(y)) %>% 
    tidyr::separate(id, c("dataset", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(n_samples = factor(n_samples, 
                                     levels = as.character(sort(unique(as.numeric(as.character(n_samples)))))))
  rownames(annotation_row) <- annotation_row$id
  
  pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", main = "Relative TPR",
           display_numbers = TRUE, colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = seq(0, 1, length.out = 101), 
           annotation_row = dplyr::select(annotation_row, n_samples, dataset), show_rownames = FALSE,
           annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
           annotation_colors = list(method = structure(cols, names = gsub(exts, "", names(cols)))[colnames(y)]),
           annotation_names_col = FALSE)
  dev.off()
  
  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_performance_relfprtpr", exts, dtpext, ".pdf"),
      width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/results_relativetruth/", ds, exts, "_results_relativetruth_summary_data.rds"))
  })
  ## Relative FPR
  y <- lapply(summary_data_list, function(m) {
    m$fpr_relative %>% 
      tidyr::separate(basemethod, c("method", "n_samples", "repl"), sep = "\\.")# %>%
  })
  y <- do.call(rbind, y) %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples))))))
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  
  ## Remove largest sample size for each data set
  for (ds in unique(y$dataset)) {
    maxc <- max(as.numeric(as.character(y$n_samples))[y$dataset == ds])
    y <- y %>% dplyr::filter(!(n_samples == maxc & dataset == ds))
  }
  
  print(ggplot(y, aes(x = method, y = FPR, color = method)) + 
          geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = n_samples)) + 
          theme_bw() + xlab("") + ylab("Relative FPR") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")))
  
  ## Relative TPR
  y <- lapply(summary_data_list, function(m) {
    m$tpr_relative %>% 
      tidyr::separate(basemethod, c("method", "n_samples", "repl"), sep = "\\.")# %>%
  })
  y <- do.call(rbind, y) %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples))))))
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  
  ## Remove largest sample size for each data set
  for (ds in unique(y$dataset)) {
    maxc <- max(as.numeric(as.character(y$n_samples))[y$dataset == ds])
    y <- y %>% dplyr::filter(!(n_samples == maxc & dataset == ds))
  }

  print(ggplot(y, aes(x = method, y = TPR, color = method)) + 
          geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = n_samples)) + 
          theme_bw() + xlab("") + ylab("Relative TPR") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")))
  
  dev.off()
}
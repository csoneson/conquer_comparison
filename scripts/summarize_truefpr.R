summarize_truefpr <- function(figdir, datasets, exts, dtpext, cols = cols,
                              singledsfigdir, cobradir, concordancedir, dschardir) {
  ## ---------------------------------- True FPR ------------------------------ ##
  pdf(paste0(figdir, "/summary_truefpr", exts, dtpext, "_1.pdf"),
      width = 10, height = 4 * length(datasets))
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0(singledsfigdir, "/truefpr/", ds, exts, "_truefpr_summary_data.rds"))
  })
  ## Heatmap of true FPRs (fraction of nominal p-values below 0.05)
  y <- lapply(summary_data_list, function(m) {
    m$fracpbelow0.05 %>% 
      tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
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
  
  pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", main = "True FPR",
           display_numbers = TRUE, colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = seq(0, 0.6, length.out = 101), 
           annotation_row = dplyr::select(annotation_row, n_samples, dataset), show_rownames = FALSE,
           annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
           annotation_colors = list(method = structure(cols, names = gsub(exts, "", names(cols)))[colnames(y)]),
           annotation_names_col = FALSE)
  dev.off()
  
  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_truefpr", exts, dtpext, "_2.pdf"),
      width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0(singledsfigdir, "/truefpr/", ds, exts, 
                   "_truefpr_summary_data.rds"))
  })
  y <- lapply(summary_data_list, function(m) {
    m$fracpbelow0.05 %>% 
      tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.")# %>%
  })
  y <- do.call(rbind, y) %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples))))))
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  
  print(ggplot(y, aes(x = method, y = FPR, color = method)) + 
          geom_hline(yintercept = 0.05) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = n_samples)) + 
          theme_bw() + xlab("") + ylab("True FPR (fraction of genes with p < 0.05)") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")))
  
  ## P-value distributions
  for (ds in datasets) {
    cbr <- readRDS(paste0(cobradir, "/", ds, exts, "_cobra.rds"))
    pv <- pval(cbr)
    tmp <- reshape2::melt(pv) %>% 
      tidyr::separate(variable, into = c("method", "ncells", "repl"), sep = "\\.") %>%
      dplyr::mutate(ncells.repl = paste0(ncells, ".", repl))
    ## Remove extension from method name
    tmp$method <- gsub(exts, "", tmp$method)
    
    for (i in unique(tmp$ncells.repl)) {
      print(tmp %>% subset(ncells.repl == i) %>% 
              ggplot(aes(x = value, fill = method)) + geom_histogram() + 
              facet_wrap(~method, scales = "free_y") + 
              theme_bw() + xlab("p-value") + ylab("") + 
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank()) + 
              scale_fill_manual(values = structure(cols, names = gsub(exts, "", names(cols))), 
                                name = "", guide = FALSE) + 
              ggtitle(paste0(ds, ".", i)))
    }
  }
  dev.off()
}
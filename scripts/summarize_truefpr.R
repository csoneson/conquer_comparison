summarize_truefpr <- function(figdir, datasets, exts, dtpext, cols,
                              singledsfigdir, cobradir, concordancedir, 
                              dschardir, origvsmockdir, distrdir, plotmethods, 
                              dstypes, pch_ncells) {
  
  gglayers <- list(
    geom_hline(yintercept = 0.05),
    geom_boxplot(outlier.size = -1),
    geom_point(position = position_jitter(width = 0.2), size = 0.5), 
    theme_bw(),
    xlab(""),
    ylab("FPR (fraction of genes with p < 0.05)"),
    scale_color_manual(values = cols),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)),
    guides(color = guide_legend(ncol = 2, title = ""))
  )
  
  ## Initialize list to hold all plots
  plots <- list()
  
  ## Read all FPR information
  truefpr <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/truefpr/", ds, e, 
                     "_truefpr_summary_data.rds"))$fracpbelow0.05
    }))
  })) %>%
    tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples)))))) %>%
    dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
    dplyr::filter(method %in% plotmethods) %>%
    dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::mutate(plot_color = cols[as.character(method)])
  
  ## Heatmap of FPRs (fraction of nominal p-values below 0.05)
  pdf(paste0(figdir, "/summary_truefpr", dtpext, "_1.pdf"),
      width = 15, height = 4 * length(datasets))
  
  for (f in unique(truefpr$filt)) {
    y <- truefpr %>% 
      dplyr::filter(filt == f) %>%
      dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", n_samples, ".", repl)) %>%
      dplyr::select(method, dataset, FPR) %>% 
      reshape2::dcast(dataset ~ method, value.var = "FPR") %>%
      tidyr::separate(dataset, c("ds", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
      dplyr::arrange(ds, as.numeric(as.character(n_samples))) %>% 
      dplyr::select(-ds, -filt, -n_samples, -repl)  %>% as.data.frame()
    rownames(y) <- y$dataset
    y$dataset <- NULL
  
    annotation_row = data.frame(id = rownames(y)) %>% 
      tidyr::separate(id, c("dataset", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
      dplyr::mutate(n_samples = factor(n_samples, 
                                       levels = as.character(sort(unique(as.numeric(as.character(n_samples)))))))
    rownames(annotation_row) <- annotation_row$id
  
    pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", 
             main = paste0("FPR, ", f), display_numbers = TRUE, 
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
             breaks = seq(0, 1, length.out = 101), 
             annotation_row = dplyr::select(annotation_row, n_samples, dataset), 
             show_rownames = FALSE, 
             annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
             annotation_colors = list(method = cols[colnames(y)]),
             annotation_names_col = FALSE)
  }
  dev.off()
  
  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_truefpr", dtpext, "_2.pdf"),
      width = 10, height = 7)

  for (f in unique(truefpr$filt)) {
    plots[[paste0("truefpr_sep_", f)]] <- 
      ggplot(truefpr %>% dplyr::filter(filt == f) %>% 
               dplyr::mutate(method = forcats::fct_reorder(method, FPR, 
                                                           fun = function(x) median(x, na.rm = TRUE),
                                                           .desc = TRUE)),
             aes(x = method, y = FPR, color = method)) + 
      gglayers + ggtitle(f)
    print(plots[[paste0("truefpr_sep_", f)]])
    
    print(plots[[paste0("truefpr_sep_", f)]] + facet_wrap(~ dtype, ncol = 1))
  }
  
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/truefpr_final", dtpext, ".pdf"), width = 12, height = 6)
  p <- plot_grid(plots$truefpr_sep_ + theme(legend.position = "none") + 
                   ggtitle("Without filtering") + scale_y_sqrt(), 
                 plots$truefpr_sep_TPM_1_25p + theme(legend.position = "none") + 
                   ggtitle("After filtering") + scale_y_sqrt(),
                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
  print(p)
  dev.off()

  pdf(paste0(figdir, "/truefpr_final", dtpext, "_bydtype.pdf"), width = 12, height = 7)
  p <- plot_grid(plots$truefpr_sep_ + facet_wrap(~ dtype, ncol = 1) + theme(legend.position = "none") + 
                   ggtitle("Without filtering") + scale_y_sqrt(), 
                 plots$truefpr_sep_TPM_1_25p + facet_wrap(~ dtype, ncol = 1) + 
                   theme(legend.position = "none") + 
                   ggtitle("After filtering") + scale_y_sqrt(),
                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
  print(p)
  dev.off()
  
  truefpr
}
summarize_truefpr <- function(figdir, datasets, exts, dtpext, cols,
                              singledsfigdir, cobradir, concordancedir, 
                              dschardir, origvsmockdir, plotmethods) {
  
  ## Generate list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_truefpr", dtpext, "_1.pdf"),
      width = 10, height = 4 * length(datasets))
  
  ## Read all true FPR information
  truefpr <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/truefpr/", ds, e, 
                     "_truefpr_summary_data.rds"))$fracpbelow0.05
    }))
  }))
  
  cols <- structure(cols, names = gsub(paste(exts, collapse = "|"), "", names(cols)))
  
  ## Heatmap of true FPRs (fraction of nominal p-values below 0.05)
  for (f in unique(truefpr$filt)) {
    y <- truefpr %>% 
      dplyr::filter(filt == f) %>%
      tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
      dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
      dplyr::filter(method %in% plotmethods) %>% 
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
             main = paste0("True FPR, ", f), display_numbers = TRUE, 
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

  truefpr <- truefpr %>% 
    tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples)))))) %>%
    dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
    dplyr::filter(method %in% plotmethods)
  
  for (f in unique(truefpr$filt)) {
    tmp <- truefpr %>% dplyr::filter(filt == f) %>%
      dplyr::group_by(method) %>% dplyr::mutate(FPRmedian = median(FPR)) %>%
      dplyr::ungroup()
    tmp$method <- factor(tmp$method, levels = unique(tmp$method[order(tmp$FPRmedian, decreasing = TRUE)]))
    plots[[paste0("truefpr_sep_", f)]] <- 
      ggplot(tmp,
             aes(x = method, y = FPR, color = method)) + 
      geom_hline(yintercept = 0.05) + geom_boxplot(outlier.size = -1) + 
      geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = n_samples)) + 
      theme_bw() + xlab("") + ylab("True FPR (fraction of genes with p < 0.05)") + 
      scale_color_manual(values = cols) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
      ggtitle(f)
    print(plots[[paste0("truefpr_sep_", f)]])
  }
  
  plots[["truefpr_comb"]] <- 
    ggplot(truefpr,
           aes(x = method, y = FPR, color = method)) + 
    geom_hline(yintercept = 0.05) + geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = n_samples)) + 
    theme_bw() + xlab("") + ylab("True FPR (fraction of genes with p < 0.05)") + 
    facet_wrap(~ filt) + 
    scale_color_manual(values = cols) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 2, title = "Number of \ncells per group"))
  print(plots[["truefpr_comb"]])
  
  ## P-value distributions
  for (ds in datasets) {
    for (e in exts) {
      cbr <- readRDS(paste0(cobradir, "/", ds, e, "_cobra.rds"))
      pv <- pval(cbr)
      tmp <- reshape2::melt(pv) %>% 
        tidyr::separate(variable, into = c("method", "ncells", "repl"), sep = "\\.") %>%
        dplyr::filter(method %in% paste0(plotmethods, e)) %>% 
        dplyr::mutate(ncells.repl = paste0(ncells, ".", repl))
      ## Remove extension from method name
      tmp$method <- gsub(paste(exts, collapse = "|"), "", tmp$method)
      
      for (i in unique(tmp$ncells.repl)) {
        p <- tmp %>% subset(ncells.repl == i) %>% 
          ggplot(aes(x = value, fill = method)) + geom_histogram() + 
          facet_wrap(~method, scales = "free_y") + 
          theme_bw() + xlab("p-value") + ylab("") + 
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) + 
          scale_fill_manual(values = structure(cols, names = gsub(paste(exts, collapse = "|"),
                                                                  "", names(cols))), 
                            name = "", guide = FALSE) + 
          ggtitle(paste0(ds, e, ".", i))
        print(p)
        if (ds == "EMTAB2805mock" & i == "48.1")
          plots[[paste0("pvalues_", ds, e, "_", i)]] <- p
      }
    }
  }
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/truefpr_final", dtpext, ".pdf"), width = 12, height = 6)
  p <- plot_grid(plot_grid(plots$truefpr_sep_ + theme(legend.position = "none") + 
                             ggtitle("Without filtering"), 
                           plots$truefpr_sep_TPM_1_25p + theme(legend.position = "none") + 
                             ggtitle("After filtering"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 get_legend(plots$truefpr_sep_ + 
                              theme(legend.position = "bottom") + 
                              guides(colour = FALSE,
                                     shape = 
                                       guide_legend(nrow = 1,
                                                    title = "Number of cells per group",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 10,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 rel_heights = c(1.7, 0.1), ncol = 1)
  print(p)
  dev.off()
  
  if ("pvalues_EMTAB2805mock_48.1" %in% names(plots)) {
    pdf(paste0(figdir, "/truefpr_final_pval_nofilt", dtpext, ".pdf"), width = 10, height = 7)
    print(plots[["pvalues_EMTAB2805mock_48.1"]] + 
            ggtitle("EMTAB2805null.48.1"))
    dev.off()
  }
  
  if ("pvalues_EMTAB2805mock_TPM_1_25p_48.1" %in% names(plots)) {
    pdf(paste0(figdir, "/truefpr_final_pval_withfilt", dtpext, ".pdf"), width = 10, height = 7)
    print(plots[["pvalues_EMTAB2805mock_TPM_1_25p_48.1"]] + 
            ggtitle("EMTAB2805null_TPM_1_25p.48.1"))
    dev.off()
  }  
  
  plots[c("truefpr_sep_", "truefpr_sep_TPM_1_25p")]
}
summarize_trueperformance <- function(figdir, datasets, exts, dtpext, cols = cols,
                                      singledsfigdir, cobradir, concordancedir, dschardir) {
  ## ---------------------------------- True FPR ------------------------------ ##
  pdf(paste0(figdir, "/summary_trueperformance", exts, dtpext, "_1.pdf"),
      width = 10, height = 4 * length(datasets))
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0(singledsfigdir, "/performance_realtruth/", ds, exts,
                   "_performance_realtruth_summary_data.rds"))
  })
  for (asp in c("FDR", "TPR")) {
    ## Heatmap of true FDRs and TPRs at padj=0.05 threshold
    y <- lapply(summary_data_list, function(m) {
      m$FDRTPR %>% 
        dplyr::filter(thr == "thr0.05") %>%
        tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
        dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", n_samples, ".", repl)) %>%
        dplyr::select_("method", "dataset", asp)
    })
    y <- do.call(rbind, y) %>% dcast(dataset ~ method, value.var = asp) %>%
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
    
    pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", 
             main = paste0("True ", asp, " at adj.p=0.05 cutoff"),
             display_numbers = TRUE, colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
             breaks = seq(0, 1, length.out = 101), 
             annotation_row = dplyr::select(annotation_row, n_samples, dataset), show_rownames = FALSE,
             annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
             annotation_colors = list(method = structure(cols, names = gsub(exts, "", 
                                                                            names(cols)))[colnames(y)]),
             annotation_names_col = FALSE)
  }
  dev.off()

  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_trueperformance", exts, dtpext, "_2.pdf"),
      width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0(singledsfigdir, "/performance_realtruth/", ds, exts,
                   "_performance_realtruth_summary_data.rds"))
  })
  y <- lapply(summary_data_list, function(m) {
    m$FDRTPR %>% 
      tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.")# %>%
  })
  y <- do.call(rbind, y) %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples))))))
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  
  for (asp in c("FDR", "TPR")) {
    p1 <- y %>% dplyr::filter(thr == "thr0.05") %>%
      ggplot(aes_string(x = "method", y = asp, color = "method")) + 
      geom_boxplot(outlier.size = -1) + 
      geom_point(position = position_jitter(width = 0.2), aes(shape = n_samples)) + 
      theme_bw() + xlab("") + ylab(paste0("True ", asp, " at adj.p = 0.05 cutoff")) + 
      scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 2, title = "Number of \ncells per group"))
    if (asp == "FDR") p1 <- p1 + geom_hline(yintercept = 0.05)
    print(p1)
    
    p2 <- y %>% dplyr::filter(thr == "thr0.05") %>%
      dplyr::mutate(ncells = paste0(n_samples, " cells per group")) %>%
      dplyr::mutate(ncells = factor(ncells, levels = paste0(sort(unique(as.numeric(as.character(gsub(" cells per group", "", ncells))))), " cells per group"))) %>%
      ggplot(aes_string(x = "ncells", y = asp, color = "method", group = "method")) + 
      geom_point(alpha = 0.25) + geom_smooth(se = FALSE) + 
      theme_bw() + xlab("") + ylab(paste0("True ", asp, " at adj.p = 0.05 cutoff")) + 
      scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 2, title = "Number of \ncells per group"))
    if (asp == "FDR") p2 <- p2 + geom_hline(yintercept = 0.05)
    print(p2)
    
    p3 <- y %>% dplyr::filter(thr == "thr0.05") %>%
      dplyr::mutate(ncells = paste0(n_samples, " cells per group")) %>%
      dplyr::mutate(ncells = factor(ncells, levels = paste0(sort(unique(as.numeric(as.character(gsub(" cells per group", "", ncells))))), " cells per group"))) %>%
      ggplot(aes_string(x = "ncells", y = asp, color = "method", group = "method")) + 
      geom_point(alpha = 0.25) + geom_smooth(se = FALSE) + 
      facet_wrap(~dataset, scales = "free_x") + 
      theme_bw() + xlab("") + ylab(paste0("True ", asp, " at adj.p = 0.05 cutoff")) + 
      scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 2, title = "Number of \ncells per group"))
    if (asp == "FDR") p3 <- p3 + geom_hline(yintercept = 0.05)
    print(p3)
  }  

  ## AUROC
  y <- lapply(summary_data_list, function(m) {
    m$AUROC %>% 
      tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.")# %>%
  })
  y <- do.call(rbind, y) %>%
    dplyr::mutate(n_samples = factor(n_samples, levels = sort(unique(as.numeric(as.character(n_samples))))))
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  
  asp <- "AUROC"
  print(y %>% 
          ggplot(aes_string(x = "method", y = asp, color = "method")) + 
          geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) + 
          theme_bw() + xlab("") + ylab("Area under ROC curve") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 1, title = "")))

  print(y %>% 
          dplyr::mutate(ncells = paste0(n_samples, " cells per group")) %>%
          dplyr::mutate(ncells = factor(ncells, levels = paste0(sort(unique(as.numeric(as.character(gsub(" cells per group", "", ncells))))), " cells per group"))) %>%
          ggplot(aes_string(x = "ncells", y = asp, color = "method", group = "method")) + 
          geom_point(alpha = 0.25) + geom_smooth(se = FALSE) + 
          theme_bw() + xlab("") + ylab("Area under ROC curve") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = "")))

  print(y %>% 
          dplyr::mutate(ncells = paste0(n_samples, " cells per group")) %>%
          dplyr::mutate(ncells = factor(ncells, levels = paste0(sort(unique(as.numeric(as.character(gsub(" cells per group", "", ncells))))), " cells per group"))) %>%
          ggplot(aes_string(x = "ncells", y = asp, color = "method", group = "method")) + 
          geom_point(alpha = 0.25) + geom_smooth(se = FALSE) + 
          facet_wrap(~dataset, scales = "free_x") + 
          theme_bw() + xlab("") + ylab("Area under ROC curve") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = "")))

  dev.off()
}
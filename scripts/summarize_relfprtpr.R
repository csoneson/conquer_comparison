summarize_relfprtpr <- function(figdir, datasets, exts, dtpext, cols,
                                singledsfigdir, cobradir, concordancedir, 
                                dschardir, origvsmockdir, distrdir, plotmethods, 
                                dstypes, pch_ncells) {

  gglayers <- list(
    geom_boxplot(outlier.size = -1),
    geom_point(position = position_jitter(width = 0.2)),
    theme_bw(),
    xlab(""),
    scale_color_manual(values = cols),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)),
    guides(color = guide_legend(ncol = 2, title = ""))
  )
  
  ## Initialize list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_relfprtpr", dtpext, "_1.pdf"),
      width = 15, height = 4 * length(datasets))
  
  X <- list()
  ## Read all relative FPR information
  X[["FPR"]] <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/results_relativetruth/", ds, e, 
                     "_results_relativetruth_summary_data.rds"))$fpr_relative
    }))
  }))
  ## Read all relative TPR information
  X[["TPR"]] <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/results_relativetruth/", ds, e, 
                     "_results_relativetruth_summary_data.rds"))$tpr_relative
    }))
  }))
  X <- lapply(X, function(x) {
    x %>%
      tidyr::separate(basemethod, c("method", "ncells", "repl"), sep = "\\.") %>%
      dplyr::mutate(ncells = factor(ncells, 
                                    levels = sort(unique(as.numeric(as.character(ncells)))))) %>%
      dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
      dplyr::filter(method %in% plotmethods)
  })  
  
  ## Heatmap of true FPRs/TPRs
  for (asp in c("FPR", "TPR")) {
    for (f in unique(X[[asp]]$filt)) {
      y <- X[[asp]] %>% 
        dplyr::filter(filt == f) %>%
        dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", ncells, ".", repl)) %>%
        dplyr::select_("method", "dataset", asp) %>% 
        reshape2::dcast(dataset ~ method, value.var = asp) %>%
        tidyr::separate(dataset, c("ds", "filt", "ncells", "repl"), sep = "\\.", remove = FALSE) %>%
        dplyr::arrange(ds, as.numeric(as.character(ncells))) %>% 
        dplyr::select(-ds, -filt, -ncells, -repl)  %>% as.data.frame()
      rownames(y) <- y$dataset
      y$dataset <- NULL
      
      annotation_row = data.frame(id = rownames(y)) %>% 
        tidyr::separate(id, c("dataset", "filt", "ncells", "repl"), sep = "\\.", remove = FALSE) %>%
        dplyr::mutate(
          ncells = factor(ncells, 
                          levels = as.character(sort(unique(as.numeric(as.character(ncells)))))))
      rownames(annotation_row) <- annotation_row$id
      
      pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", 
               main = paste0("Relative ", asp, ", ", f), display_numbers = TRUE, 
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               breaks = seq(0, 1, length.out = 101), 
               annotation_row = dplyr::select(annotation_row, ncells, dataset), 
               show_rownames = FALSE, 
               annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
               annotation_colors = list(method = cols[colnames(y)]),
               annotation_names_col = FALSE)
    }
  }  
  dev.off()
  
  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_relfprtpr", dtpext, "_2.pdf"),
      width = 10, height = 7)
  
  X <- lapply(X, function(x) {
    ## Remove largest sample size for each data set
    for (ds in unique(x$dataset)) {
      for (f in unique(x$filt)) {
        maxc <- max(as.numeric(as.character(x$ncells))[x$dataset == ds & x$filt == f])
        x <- x %>% dplyr::filter(!(ncells == maxc & dataset == ds & filt == f))
      }
    }
    x
  })  
  
  for (asp in c("FPR", "TPR")) {
    for (f in unique(X[[asp]]$filt)) {
      plots[[paste0(asp, "_sep_", f)]] <- 
        ggplot(X[[asp]] %>% dplyr::filter(filt == f),
               aes_string(x = "method", y = asp, color = "method")) + 
        gglayers + ylab(paste0("Relative ", asp)) + 
        ggtitle(f)
      print(plots[[paste0(asp, "_sep_", f)]])
    }
  }  
  dev.off()
}
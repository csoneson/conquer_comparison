summarize_nbrdet <- function(figdir, datasets, exts, dtpext, cols,
                             singledsfigdir, cobradir, concordancedir, 
                             dschardir, origvsmockdir) {
  
  ## Generate list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_nbrdet", dtpext, ".pdf"), width = 14, height = 7)
  
  ## Read all necessary information, for all filterings
  nbrgenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, 
                     "_nbr_called.rds")) %>%
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method))
    }))
  }))
  nbrgenes <- nbrgenes %>% 
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells))))
  
  ## Set plot symbols for number of cells per group
  ncells <- levels(nbrgenes$ncells_fact)
  pch <- c(16, 17, 15, 3, 7, 8, 4, 6, 9, 10, 11, 12, 13, 14)[1:length(ncells)]
  names(pch) <- as.character(ncells)
  ## Define colors for plotting
  cols <- structure(cols, names = gsub(paste(exts, collapse = "|"), "", names(cols)))
  
  for (f in unique(nbrgenes$filt)) {
    plots[[paste0("nbrdet_sep_", f)]] <- 
      ggplot(nbrgenes %>% dplyr::filter(filt == f), 
             aes(x = method, y = nbr_sign0.05, color = method)) + 
      geom_boxplot(outlier.size = -1) +
      geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
      theme_bw() + xlab("") + ylab("Number of genes with adjusted p-value below 0.05") + 
      scale_color_manual(values = cols) + 
      scale_shape_manual(values = pch) + 
      facet_wrap(~dataset, scales = "free_y") + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      ggtitle(f)
    print(plots[[paste0("nbrdet_sep_", f)]])
    
    ## Line plot
    plots[[paste0("nbrdet_sep_line_", f)]] <- 
      ggplot(nbrgenes %>% dplyr::filter(filt == f), 
             aes(x = ncells_fact, y = nbr_sign0.05, group = method, color = method)) + 
      geom_point(alpha = 0.25, size = 1) + geom_smooth(size = 0.75, se = FALSE, 
                                                       method = "loess", span = 1) + 
      theme_bw() + xlab("") + ylab("Number of genes with adjusted p-value below 0.05") + 
      scale_color_manual(values = cols) + 
      facet_wrap(~dataset, scales = "free") + 
      guides(color = guide_legend(ncol = 2, title = "")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      ggtitle(f)
    print(plots[[paste0("nbrdet_sep_line_", f)]])
  }
  
  plots[["nbrdet_comb"]] <- 
    ggplot(nbrgenes %>% dplyr::group_by(dataset, filt, ncells_fact, repl) %>%
             dplyr::mutate(nbr_sign0.05_rel = nbr_sign0.05/max(nbr_sign0.05)), 
           aes(x = method, y = nbr_sign0.05_rel, color = method)) + 
    geom_boxplot(outlier.size = -1) +
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
    facet_wrap(~ filt, nrow = 1) + 
    theme_bw() + xlab("") + ylab("Relative number of genes \nwith adjusted p-value below 0.05") + 
    scale_color_manual(values = cols) + 
    scale_shape_manual(values = pch) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  print(plots[["nbrdet_comb"]])
  
  for (f in unique(nbrgenes$filt)) {
    plots[[paste0("nbrdet_comb_", f)]] <-
      ggplot(nbrgenes %>% dplyr::filter(filt == f) %>% 
               dplyr::group_by(dataset, filt, ncells_fact, repl) %>%
               dplyr::mutate(nbr_sign0.05_rel = nbr_sign0.05/max(nbr_sign0.05)), 
             aes(x = method, y = nbr_sign0.05_rel, color = method)) + 
      geom_boxplot(outlier.size = -1) +
      geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
      theme_bw() + xlab("") + ylab("Relative number of genes \nwith adjusted p-value below 0.05") + 
      scale_color_manual(values = cols) + 
      scale_shape_manual(values = pch) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      ggtitle(f)
    print(plots[[paste0("nbrdet_comb_", f)]])
    
    tmp <- nbrgenes %>% dplyr::filter(filt == f) %>% 
      dplyr::group_by(dataset, filt, ncells_fact, repl) %>%
      dplyr::mutate(nbr_sign0.05_rel = nbr_sign0.05/max(nbr_sign0.05)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(method) %>%
      dplyr::mutate(nbr_sign0.05_rel_median = median(nbr_sign0.05_rel)) %>%
      dplyr::ungroup()
    tmp$method <- factor(tmp$method, levels = unique(tmp$method[order(tmp$nbr_sign0.05_rel_median, 
                                                                      decreasing = TRUE)]))
    plots[[paste0("nbrdet_comb_", f, "_sorted")]] <-
      ggplot(tmp, 
             aes(x = method, y = nbr_sign0.05_rel, color = method)) + 
      geom_boxplot(outlier.size = -1) +
      geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
      theme_bw() + xlab("") + ylab("Relative number of genes \nwith adjusted p-value below 0.05") + 
      scale_color_manual(values = cols) + 
      scale_shape_manual(values = pch) + 
      guides(color = guide_legend(ncol = 2, title = ""),
             shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      ggtitle(f)
    print(plots[[paste0("nbrdet_comb_", f, "_sorted")]])
  }
  
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/nbrdet_final", dtpext, ".pdf"), width = 12, height = 6)
  p <- plot_grid(plot_grid(plots[["nbrdet_comb__sorted"]] + theme(legend.position = "none") + 
                             ggtitle("Without filtering"), 
                           plots[["nbrdet_comb_TPM_1_25p_sorted"]] + theme(legend.position = "none") + 
                             ggtitle("After filtering"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 get_legend(plots[["nbrdet_comb__sorted"]] + 
                              theme(legend.position = "bottom") + 
                              guides(colour = FALSE,
                                     shape = 
                                       guide_legend(nrow = 1,
                                                    title = "Number of cells per group",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 rel_heights = c(1.7, 0.1), ncol = 1)
  print(p)
  dev.off()
  
  ## Split by data set
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_byds.pdf"), width = 18, height = 10)
  p <- plots[["nbrdet_sep_"]] + 
    theme(legend.position = "bottom") + 
    guides(colour = FALSE,
           shape = 
             guide_legend(nrow = 1,
                          title = "Number of cells per group",
                          override.aes = list(size = 1.5),
                          title.theme = element_text(size = 12, angle = 0),
                          label.theme = element_text(size = 12, angle = 0),
                          keywidth = 1, default.unit = "cm")) + 
    ggtitle("Without filtering")
  print(p)
  dev.off()
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_byds_filtered.pdf"), width = 18, height = 10)
  p <- plots[["nbrdet_sep_TPM_1_25p"]] + 
    theme(legend.position = "bottom") + 
    guides(colour = FALSE,
           shape = 
             guide_legend(nrow = 1,
                          title = "Number of cells per group",
                          override.aes = list(size = 1.5),
                          title.theme = element_text(size = 12, angle = 0),
                          label.theme = element_text(size = 12, angle = 0),
                          keywidth = 1, default.unit = "cm")) + 
    ggtitle("After filtering")
  print(p)
  dev.off()
  
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_line_byds.pdf"), width = 18, height = 10)
  p <- plots[["nbrdet_sep_line_"]] + 
    theme(legend.position = "bottom") + 
    guides(shape = FALSE,
           colour = 
             guide_legend(nrow = 3,
                          title = "",
                          override.aes = list(size = 1.5),
                          title.theme = element_text(size = 12, angle = 0),
                          label.theme = element_text(size = 12, angle = 0),
                          keywidth = 1, default.unit = "cm")) + 
    ggtitle("Without filtering")
  print(p)
  dev.off()
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_line_byds_filtered.pdf"), width = 18, height = 10)
  p <- plots[["nbrdet_sep_line_TPM_1_25p"]] + 
    theme(legend.position = "bottom") + 
    guides(shape = FALSE,
           colour = 
             guide_legend(nrow = 3,
                          title = "",
                          override.aes = list(size = 1.5),
                          title.theme = element_text(size = 12, angle = 0),
                          label.theme = element_text(size = 12, angle = 0),
                          keywidth = 1, default.unit = "cm")) + 
    ggtitle("After filtering")
  print(p)
  dev.off()
}
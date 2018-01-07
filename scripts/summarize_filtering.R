summarize_filtering <- function(figdir, datasets, exts, dtpext, cols,
                                singledsfigdir, cobradir, concordancedir, 
                                dschardir, origvsmockdir, distrdir, plotmethods, 
                                dstypes, pch_ncells) {
  exts <- setdiff(exts, "")
  
  for (e in exts) {
    pdf(paste0(figdir, "/summary_filtering", e, dtpext, ".pdf"), width = 12, height = 7)
    
    ## Read information about number of tested genes before and after filtering
    summary_data_list_orig <- lapply(datasets, function(ds) {
      readRDS(paste0(cobradir, "/", ds, "_nbr_called.rds"))
    })
    summary_data_list_filt <- lapply(datasets, function(ds) {
      readRDS(paste0(cobradir, "/", ds, e, "_nbr_called.rds"))
    })
    
    ## Merge data across data sets
    L <- rbind(do.call(rbind, summary_data_list_orig),
               do.call(rbind, summary_data_list_filt)) %>%
      dplyr::mutate(method = gsub(e, "", method)) %>%
      dplyr::filter(method == method[1]) %>%
      dplyr::group_by(dataset, ncells, repl) %>%
      dplyr::summarize(filtered = nbr_tested[filt == gsub("^_", "", e)],
                       unfiltered = nbr_tested[filt == ""],
                       retain = filtered/unfiltered) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(ncells = factor(ncells, levels = sort(unique(as.numeric(as.character(ncells))))))
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    ncells_col <- gg_color_hue(n = length(unique(L$ncells)))
    names(ncells_col) <- sort(unique(L$ncells))
    
    p1 <- ggplot(L, aes(x = dataset, y = retain)) + geom_boxplot(outlier.size = -1) + 
      geom_point(position = position_jitter(width = 0.2), aes(color = ncells)) + 
      theme_bw() + xlab("") + 
      scale_color_manual(values = ncells_col) + 
      ylab("Retained fraction of original\nset of genes after filtering") + 
      guides(color = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      stat_summary(fun.data = function(x) {
        return(data.frame(y = 0.66,
                          label = paste0("n=", sum(!is.na(x)))))}, 
        geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
        hjust = 1, angle = 90) + 
      geom_hline(yintercept = 0.6, linetype = "dashed")
    
    p2 <- ggplot(L %>% dplyr::select(-retain) %>%
                   tidyr::gather(dst, nbrgenes, filtered, unfiltered) %>% 
                   dplyr::mutate(dst = factor(dst, levels = c("unfiltered", "filtered"))),
                 aes(x = dataset, y = nbrgenes, fill = dst)) + 
      geom_boxplot(outlier.size = -1, alpha = 0.3) + 
      geom_point(aes(color = ncells, group = dst), 
                 position = position_jitterdodge(jitter.width = 0.7)) + 
      scale_color_manual(values = ncells_col) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13)) + 
      ylab("Number of retained genes") + xlab("") + 
      guides(color = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
      scale_fill_manual(values = c(filtered = "lightblue", unfiltered = "pink")) + 
      geom_vline(xintercept = 0.5 + seq_len(length(unique(L$dataset)) - 1),
                 linetype = "dashed", alpha = 0.3) + 
      guides(alpha = FALSE, fill = FALSE) + 
      stat_summary(fun.data = function(x) {
        return(data.frame(y = 55000,
                          label = paste0("n=", sum(!is.na(x)))))}, 
        geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
        hjust = 1, angle = 90) + 
      geom_hline(yintercept = 50000, linetype = "dashed")
    
    pcomb <- plot_grid(plot_grid(p1 + theme(legend.position = "none"), 
                                 p2 + theme(legend.position = "none"),
                                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                       get_legend(p1 + theme(legend.position = "bottom") + 
                                    guides(color = 
                                             guide_legend(nrow = 2,
                                                          title = "Number of cells per group",
                                                          override.aes = list(size = 1.5),
                                                          title.theme = element_text(size = 12,
                                                                                     angle = 0),
                                                          label.theme = element_text(size = 10,
                                                                                     angle = 0),
                                                          keywidth = 1, default.unit = "cm"))),
                       rel_heights = c(1.7, 0.15), ncol = 1)
    print(pcomb)
    
    ## Plot library size vs fraction zeros across all data sets
    summary_data_list <- lapply(datasets, function(ds) {
      readRDS(paste0(dschardir, "/", ds, "_dataset_characteristics_summary_data.rds"))
    })
    L <- do.call(
      rbind, 
      lapply(summary_data_list, function(x) x$char_cells_m %>% 
               dplyr::filter(mtype %in% c("libsize", "fraczeroround")) %>%
               dplyr::mutate(cell = paste0(dataset, ".", cell, ".", ncells, ".", repl)) %>%
               dplyr::select(cell, mtype, value) %>%
               tidyr::spread(key = mtype, value = value) %>%
               tidyr::separate(cell, into = c("dataset", "cell", "ncells", "repl"), sep = "\\."))) %>%
      dplyr::mutate(
        ncells = factor(ncells, 
                        levels = paste0(sort(unique(
                          as.numeric(as.character(gsub(" cells per group", "", ncells))))), 
                          " cells per group")))
    pp <- ggplot(L, aes(x = libsize, y = fraczeroround, color = dataset, shape = ncells)) + 
      geom_point() + theme_bw() + xlab("Library size") + 
      ylab("Fraction zeros after rounding") + 
      scale_shape_manual(values = seq_len(length(unique(L$ncells))), name = "") + 
      scale_color_discrete(name = "") + 
      theme(axis.text = element_text(size = 12), 
            axis.title = element_text(size = 13))
    print(pp)
    print(pp + scale_x_log10())

    dev.off()
    
    pdf(paste0(figdir, "/filtering_final", e, dtpext, ".pdf"), width = 12, height = 7)
    print(pcomb)
    dev.off()
    
    if (dtpext == "_real") {
      pdf(paste0(figdir, "/filtering_for_slides", e, dtpext, ".pdf"), width = 10, height = 5.6)
      print(plot_grid(p1 + theme(legend.position = "none"), 
                      p2 + theme(legend.position = "none"),
                      align = "h", rel_widths = c(1, 1), nrow = 1))
      dev.off()
    }
  }
}
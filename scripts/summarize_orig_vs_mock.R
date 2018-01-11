summarize_orig_vs_mock <- function(figdir, datasets, exts, dtpext, cols,
                                   singledsfigdir, cobradir, concordancedir, 
                                   dschardir, origvsmockdir, distrdir, plotmethods, 
                                   dstypes, pch_ncells) {

  gglayers <- list(
    geom_boxplot(outlier.size = -1),
    theme_bw(),
    xlab(""),
    scale_color_manual(values = structure(cols, names = gsub(paste(exts, collapse = "|"),
                                                             "", names(cols))), name = ""),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          legend.position = "none"),
    geom_point(position = position_jitter(width = 0.2), size = 0.5),
    ylim(0, 1),
    facet_wrap(~ dataset)
  )
  
  ## Initialize list to hold all plots
  plots <- list()

  pdf(paste0(figdir, "/summary_orig_vs_mock", dtpext, ".pdf"), 
      width = 10, height = 7)
  
  concordances <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(origvsmockdir, "/", ds, e, 
                     "_orig_vs_mock_summary_data.rds"))$concordances %>%
        dplyr::ungroup() %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
    }))
  })) %>% dplyr::mutate(plot_color = cols[as.character(method)]) %>%
    dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::select(method, dataset, dtype, filt, ncells, replicate1, replicate2, 
                  k, AUCs, tp)
  
  dslist <- lapply(list(all = unique(concordances$dataset),
                        sub = c("GSE60749-GPL13112", "10XMonoCytoT")), 
                   function(x) intersect(x, concordances$dataset)
  )
  (dslist <- dslist[sapply(dslist, length) > 0])
  
  for (f in unique(concordances$filt)) {
    for (k0 in unique(concordances$k)) {
      concs0 <- concordances %>% dplyr::filter(filt == f) %>%
        dplyr::filter(k == k0)
      for (nm in names(dslist)) {
        concs <- concs0  %>%
          dplyr::filter(dataset %in% dslist[[nm]]) %>%
          dplyr::group_by(method, dataset, filt, ncells, k) %>%
          dplyr::mutate(tokeep = length(unique(tp)) == 2) %>%
          dplyr::filter(tokeep) %>% 
          dplyr::mutate(signal_vs_mock = mean(AUCs[tp == "signal"]) - mean(AUCs[tp == "mock"])) %>%
          dplyr::ungroup()
        
        xtmp <- concs %>% 
          dplyr::mutate(method = forcats::fct_reorder(method, signal_vs_mock, 
                                                      fun = median, na.rm = TRUE,
                                                      .desc = TRUE))
        x0 <- xtmp %>% dplyr::filter(tp == "signal")
        p <- x0 %>%
          ggplot(aes(x = method, y = AUCs, col = method)) + 
          ylab("Area under concordance curve, signal data set") + 
          gglayers + ggtitle(paste0(f, ", top-", k0, " genes")) + 
          stat_summary(fun.data = function(x) {
            return(data.frame(y = 1.05,
                              label = paste0("n=", sum(!is.na(x)))))}, 
            geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
            hjust = -0.2, angle = 90) + 
          scale_y_continuous(limits = c(0, 1.3), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
          geom_hline(yintercept = 1.05, linetype = "dashed")
        print(p)
        plots[[paste0("auc_signal_sep_", nm, "_", f, "_", k0)]] <- p
        
        x0 <- xtmp %>% dplyr::filter(tp == "mock")
        p <- x0 %>%
          ggplot(aes(x = method, y = AUCs, col = method)) + 
          ylab("Area under concordance curve, null data set") + 
          gglayers + ggtitle(paste0(f, ", top-", k0, " genes")) + 
          stat_summary(fun.data = function(x) {
            return(data.frame(y = 1.05,
                              label = paste0("n=", sum(!is.na(x)))))}, 
            geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
            hjust = -0.2, angle = 90) + 
          scale_y_continuous(limits = c(0, 1.3), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
          geom_hline(yintercept = 1.05, linetype = "dashed")
        print(p)
        plots[[paste0("auc_mock_sep_", nm, "_", f, "_", k0)]] <- p
      }
    }
  }
  dev.off()
    
  ## -------------------------- Final summary plots ------------------------- ##

  for (k0 in unique(concordances$k)) {
    for (nm in names(dslist)) {
      ## Split by data set
      pdf(paste0(figdir, "/orig_vs_mock_final", dtpext, "_", k0, "_sepbyds_", nm, ".pdf"), 
          width = 15, height = length(dslist[[nm]]) * 2 + 3)
      print(plot_grid(ggdraw() + 
                        draw_label(paste0("After filtering, top-", k0, " genes"),
                                   fontface = "bold"), 
                      plot_grid(
                        plots[[paste0("auc_signal_sep_", nm, "_TPM_1_25p_", k0)]] + 
                          facet_wrap(~dataset, ncol = 1) + 
                          ggtitle(""),
                        plots[[paste0("auc_mock_sep_", nm, "_TPM_1_25p_", k0)]] + 
                          facet_wrap(~dataset, ncol = 1) + 
                          ggtitle(""),
                        ncol = 2, rel_widths = c(1, 1)
                      ), ncol = 1, rel_heights = c(0.25, 6)
      ))
      dev.off()
    }
  }
  
  concordances %>%
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells))))
}

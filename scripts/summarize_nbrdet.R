summarize_nbrdet <- function(figdir, datasets, exts, dtpext, cols,
                             singledsfigdir, cobradir, concordancedir, 
                             dschardir, origvsmockdir, distrdir, plotmethods, 
                             dstypes, pch_ncells) {
  
  gglayers <- list(
    theme_bw(),
    ylab("Number of genes with adjusted p-value below 0.05"),
    scale_color_manual(values = cols),
    guides(color = FALSE),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          strip.text = element_text(size = 12))
  )
  gglayersp <- c(list(geom_boxplot(outlier.size = -1),
                      geom_point(position = position_jitter(width = 0.2), size = 0.5),
                      xlab("")),
                 gglayers)
  gglayersl <- c(list(geom_point(alpha = 0.25, size = 1),
                      geom_smooth(size = 0.75, se = FALSE, method = "loess", span = 1),
                      xlab("Number of cells per group")),
                 gglayers)
  
  ## Initialize list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_nbrdet", dtpext, ".pdf"), width = 14, height = 7)
  
  ## Read all necessary information, for all filterings
  nbrgenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, 
                     "_nbr_called.rds")) %>%
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
    }))
  })) %>% 
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells)))) %>%
    dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::mutate(plot_color = cols[as.character(method)]) %>%
    dplyr::group_by(dataset, filt, ncells_fact, repl) %>%
    dplyr::mutate(nbr_sign_adjp0.05_rel = nbr_sign_adjp0.05/max(nbr_sign_adjp0.05)) %>%
    dplyr::ungroup()
  
  ## ------------------------------------------------------------------------ ##
  
  for (f in unique(nbrgenes$filt)) {
    x0 <- nbrgenes %>% dplyr::filter(filt == f) %>% 
      dplyr::mutate(method = forcats::fct_reorder(method, nbr_sign_adjp0.05_rel, 
                                                  fun = median, na.rm = TRUE,
                                                  .desc = TRUE))
    plots[[paste0("nbrdet_sep_", f)]] <- 
      ggplot(x0, aes(x = method, y = nbr_sign_adjp0.05, color = method)) + 
      gglayersp + facet_wrap(~ dataset, scales = "fixed") + ggtitle(f) + 
      scale_y_sqrt() + 
      stat_summary(fun.data = function(x) {
        return(data.frame(y = max(x) + sqrt(1000),
                          label = paste0("n=", sum(!is.na(x)))))}, 
        geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5,
        hjust = 1, angle = 90)
    print(plots[[paste0("nbrdet_sep_", f)]])
    
    ## Line plot
    plots[[paste0("nbrdet_sep_line_", f)]] <- 
      ggplot(nbrgenes %>% dplyr::filter(filt == f), 
             aes(x = ncells_fact, y = nbr_sign_adjp0.05, group = method, color = method)) + 
      gglayersl + facet_wrap(~dataset, scales = "free") + ggtitle(f)
    print(plots[[paste0("nbrdet_sep_line_", f)]])
  }
  
  ## Compare before and after filtering
  x0 <- nbrgenes %>%
    dplyr::group_by(method, ncells_fact, repl, dataset, dtype, plot_color) %>%
    dplyr::filter(length(nbr_sign_adjp0.05) == 2) %>%
    dplyr::summarize(nbr_sign_adjp0.05_unfilt = nbr_sign_adjp0.05[filt == ""],
                     nbr_sign_adjp0.05_filt = nbr_sign_adjp0.05[filt == "TPM_1_25p"],
                     ratio_nbr_sign_adjp0.05 = 
                       nbr_sign_adjp0.05[filt == "TPM_1_25p"]/nbr_sign_adjp0.05[filt == ""]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ratio_nbr_sign_adjp0.05 = replace(ratio_nbr_sign_adjp0.05, 
                                                    ratio_nbr_sign_adjp0.05 > 10, 10)) %>%
    dplyr::mutate(method = forcats::fct_reorder(method, ratio_nbr_sign_adjp0.05, 
                                                fun = median, na.rm = TRUE,
                                                .desc = TRUE))
  plots[["nbrdet_bydtype_ratio"]] <- 
    ggplot(x0, aes(x = method, y = ratio_nbr_sign_adjp0.05, color = method)) + 
    geom_hline(yintercept = 1) + gglayersp + facet_wrap(~ dtype, ncol = 1) + scale_y_sqrt() + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = sqrt(14),
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5,
      hjust = 1, angle = 90) + 
    ylab("Ratio between number of genes with adjusted p-value\nbelow 0.05 for filtered and unfiltered data set instances")
  print(plots[["nbrdet_bydtype_ratio"]])
  
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  ## Box plots
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_byds.pdf"), width = 19, height = 10)
  print(plots[["nbrdet_sep_"]] +  
          guides(colour = FALSE) + ggtitle("Without filtering"))
  dev.off()
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_byds_filtered.pdf"), width = 19, height = 10)
  print(plots[["nbrdet_sep_TPM_1_25p"]] + 
          guides(colour = FALSE) + ggtitle("After filtering"))
  dev.off()
  
  ## Ratios
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_ratio_bydtype.pdf"), width = 9, height = 6.5)
  print(plots[["nbrdet_bydtype_ratio"]] + guides(colour = FALSE))
  dev.off()
  
  ## Line plots
  pdf(paste0(figdir, "/nbrdet_final", dtpext, "_line_byds.pdf"), width = 18, height = 10)
  p <- plots[["nbrdet_sep_line_"]] + 
    theme(legend.position = "bottom") + 
    guides(colour = 
             guide_legend(nrow = 4,
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
    guides(colour = 
             guide_legend(nrow = 4,
                          title = "",
                          override.aes = list(size = 1.5),
                          title.theme = element_text(size = 12, angle = 0),
                          label.theme = element_text(size = 12, angle = 0),
                          keywidth = 1, default.unit = "cm")) + 
    ggtitle("After filtering")
  print(p)
  dev.off()
  
  nbrgenes %>% dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, 
                             nbr_sign_adjp0.05, nbr_sign_adjp0.05_rel)
}
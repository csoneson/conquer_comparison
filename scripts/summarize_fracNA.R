summarize_fracNA <- function(figdir, datasets, exts, dtpext, cols,
                             singledsfigdir, cobradir, concordancedir, 
                             dschardir, origvsmockdir, distrdir, plotmethods, 
                             dstypes, pch_ncells) {
  
  ## Define layers to reuse across ggplots
  gglayers <- list(
    geom_boxplot(outlier.size = -1),
    geom_point(position = position_jitter(width = 0.2), size = 0.5),
    theme_bw(),
    xlab(""),
    scale_color_manual(values = cols),
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.2)),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          strip.text = element_text(size = 12),
          legend.position = "none")
  )
  
  ## Initialize list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_fracNA", dtpext, ".pdf"), width = 14, height = 7)
  
  ## Read all necessary information, for all filterings
  nbrgenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, "_nbr_called.rds")) %>%
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
    }))
  })) %>% 
    dplyr::mutate(fracNA = nbr_NA/nbr_tested) %>%
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells)))) %>%
    dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::mutate(dataset = gsub("mock", "null", dataset)) %>%
    dplyr::mutate(plot_color = cols[as.character(method)])
  
  for (f in unique(nbrgenes$filt)) {
    ## All data sets in one plot
    plots[[paste0("fracna_comb_", f)]] <-
      ggplot(nbrgenes %>% dplyr::filter(filt == f) %>% 
               dplyr::mutate(method = forcats::fct_reorder(method, fracNA, 
                                                           fun = median, na.rm = TRUE,
                                                           .desc = TRUE)),
             aes(x = method, y = fracNA, color = method)) + 
      gglayers + ylab("Fraction of NA adjusted p-values") + 
      ggtitle(f) + 
      stat_summary(fun.data = function(x) {
        return(data.frame(y = 1.05,
                          label = paste0("n=", sum(!is.na(x)))))}, 
        geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5,
        hjust = -0.2, angle = 90) + 
      geom_hline(yintercept = 1.05, linetype = "dashed")
    print(plots[[paste0("fracna_comb_", f)]])
    
    ## Separate plot for each data set
    print(plots[[paste0("fracna_comb_", f)]] + facet_wrap(~ dataset))
    
    ## Facet by dataset type
    print(plots[[paste0("fracna_comb_", f)]] + facet_wrap(~ dtype, ncol = 1))
  }
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/fracNA_final", dtpext, ".pdf"), width = 12, height = 6)
  p <- plot_grid(plots[["fracna_comb_"]] + ggtitle("Without filtering"), 
                 plots[["fracna_comb_TPM_1_25p"]] + ggtitle("After filtering"),
                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
  print(p)
  dev.off()
  
  ## Split by data type
  pdf(paste0(figdir, "/fracNA_final", dtpext, "_bydtype.pdf"), width = 12, height = 7)
  p <- plot_grid(plots[["fracna_comb_"]] + facet_wrap(~ dtype, ncol = 1) + 
                   ggtitle("Without filtering"), 
                 plots[["fracna_comb_TPM_1_25p"]] + facet_wrap(~ dtype, ncol = 1) + 
                   ggtitle("After filtering"),
                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
  print(p)
  dev.off()
  
  ## Split by data set
  pdf(paste0(figdir, "/fracNA_final", dtpext, "_byds.pdf"), width = 21, height = 22)
  p <- plots[["fracna_comb_"]] + facet_wrap(~ dataset, ncol = 4) + 
    ggtitle("Without filtering")
  print(p)
  dev.off()
  
  if (dtpext == "_real")
    write.table(nbrgenes %>% dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, fracNA) %>%
                  dplyr::mutate(dataset = gsub("mock", "null", dataset)),
                file = "export_results/Figure1.csv", row.names = FALSE,
                col.names = TRUE, sep = ",", quote = FALSE)
  
  nbrgenes %>% dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, fracNA)
}
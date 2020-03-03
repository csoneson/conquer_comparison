summarize_pvalhist <- function(figdir, datasets, exts, dtpext, cols,
                               singledsfigdir, cobradir, concordancedir, 
                               dschardir, origvsmockdir, distrdir, plotmethods, 
                               dstypes, pch_ncells) {

  ## Initalize list to hold all plots
  plots <- list()
  
  ## P-value distributions
  for (ds in c("EMTAB2805mock")) {
    for (e in exts) {
      cbr <- readRDS(paste0(cobradir, "/", ds, e, "_cobra.rds"))
      pv <- pval(cbr)
      tmp <- reshape2::melt(pv) %>%
        tidyr::separate(variable, into = c("method", "ncells", "repl"), sep = "\\.") %>%
        dplyr::filter(method %in% paste0(plotmethods, e)) %>%
        dplyr::mutate(ncells.repl = paste0(ncells, ".", repl)) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method))
      
      for (i in c("48.1")) {
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
        plots[[paste0("pvalues_", ds, e, "_", i)]] <- p
      }
    }
  }
  
  if ("pvalues_EMTAB2805mock_48.1" %in% names(plots)) {
    pdf(paste0(figdir, "/pvalhist_final_nofilt", dtpext, ".pdf"), width = 10, height = 7)
    print(plots[["pvalues_EMTAB2805mock_48.1"]] + 
            ggtitle("EMTAB2805null.48.1"))
    dev.off()
  }
  
  if ("pvalues_EMTAB2805mock_TPM_1_25p_48.1" %in% names(plots)) {
    pdf(paste0(figdir, "/pvalhist_final_filt", dtpext, ".pdf"), width = 10, height = 7)
    print(plots[["pvalues_EMTAB2805mock_TPM_1_25p_48.1"]] + 
            ggtitle("EMTAB2805null_TPM_1_25p.48.1"))
    dev.off()
  }  
  
}
  
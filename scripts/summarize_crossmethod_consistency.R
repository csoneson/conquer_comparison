source("scripts/help_function_crossmethod_concordance.R")

summarize_crossmethod_consistency <- function(figdir, datasets, exts, dtpext, cols,
                                              singledsfigdir, cobradir, concordancedir, 
                                              dschardir, origvsmockdir) {
  plots <- list()
  
  K0 <- c(100, 1000)
  
  pdf(paste0(figdir, "/summary_crossmethod_consistency", exts, dtpext, ".pdf"), 
      width = 14, height = 10)
  
  concordances <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      as.data.frame(readRDS(paste0(concordancedir, "/", ds, 
                                   e, "_concordances.rds"))$concordance_betweenmethods_pairwise) %>%
        dplyr::filter(k %in% K0) %>% 
        dplyr::mutate(method1 = gsub(paste(exts, collapse = "|"), "", method1)) %>%
        dplyr::mutate(method2 = gsub(paste(exts, collapse = "|"), "", method2))
    }))
  }))

  for (f in unique(concordances$filt)) {
    for (k0 in K0) {
      plots[[paste0(f, "_", k0)]] <- 
        help_function_crossmethod_concordance(concordances %>% dplyr::filter(filt == f), 
                                              k0 = k0, titleext = f)
    }
  }
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/crossmethod_consistency_final", dtpext, ".pdf"), 
      width = 14, height = 10)
  print(plots[["TPM_1_25p_100"]]$concordancedistr_color + 
          ggtitle("After filtering"))
  dev.off()
  
}
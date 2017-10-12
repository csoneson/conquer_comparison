source("scripts/help_function_crossmethod_concordance.R")

summarize_crossmethod_consistency <- function(figdir, datasets, exts, dtpext, cols,
                                              singledsfigdir, cobradir, concordancedir, 
                                              dschardir, origvsmockdir, distrdir,
                                              plotmethods, dstypes, pch_ncells) {
  plots <- list()
  
  K0 <- c(100, 1000)
  
  pdf(paste0(figdir, "/summary_crossmethod_consistency", dtpext, ".pdf"), 
      width = 14, height = 10)
  
  concordances <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      as.data.frame(readRDS(paste0(concordancedir, "/", ds, 
                                   e, "_concordances.rds"))$concordance_betweenmethods_pairwise) %>%
        dplyr::filter(k %in% K0) %>% 
        dplyr::mutate(method1 = gsub(paste(exts, collapse = "|"), "", method1)) %>%
        dplyr::mutate(method2 = gsub(paste(exts, collapse = "|"), "", method2)) %>%
        dplyr::filter(method1 %in% plotmethods & method2 %in% plotmethods) %>% 
        dplyr::mutate(repl = as.numeric(as.character(repl))) %>%
        dplyr::mutate(ncells = as.numeric(as.character(ncells)))
    }))
  }))
  
  dsinfo <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      as.data.frame(readRDS(paste0(dschardir, "/", ds, 
                                   e, "_dataset_characteristics_summary_data.rds"))$char_ds_m) %>%
        dplyr::rename(ncells = n_cells) %>%
        dplyr::mutate(ncells = as.numeric(as.character(ncells)))
    }))
  }))

  concordances <- dplyr::full_join(concordances, dsinfo, 
                                   by = c("ncells", "repl", "dataset", "filt"))
  
  for (f in unique(concordances$filt)) {
    for (k0 in K0) {
      plots[[paste0(f, "_", k0)]] <- 
        help_function_crossmethod_concordance(concordances %>% dplyr::filter(filt == f), 
                                              k0 = k0, titleext = f)
    }
  }
  dev.off()
  
  print(names(plots))
  
  ## -------------------------- Final summary plots ------------------------- ##
  for (k0 in K0) {
    pdf(paste0(figdir, "/crossmethod_consistency_final", dtpext, "_", k0, ".pdf"), 
        width = 12, height = 12)
    print(plots[[paste0("TPM_1_25p_", k0)]]$concordancedistr_color + 
            ggtitle("After filtering"))
    dev.off()
    
    if ("concordance_dep_ncells_color" %in% names(plots[[paste0("TPM_1_25p_", k0)]])) {
      pdf(paste0(figdir, "/crossmethod_consistency_final", dtpext, "_", k0, "_ncellsdep.pdf"), 
          width = 12, height = 12)
      print(plots[[paste0("TPM_1_25p_", k0)]]$concordance_dep_ncells_color + 
              ggtitle("After filtering"))
      dev.off()
    }
    
    if ("concordance_dep_silhouette_color" %in% names(plots[[paste0("TPM_1_25p_", k0)]])) {
      pdf(paste0(figdir, "/crossmethod_consistency_final", dtpext, "_", k0, "_silhouettedep.pdf"), 
          width = 12, height = 12)
      print(plots[[paste0("TPM_1_25p_", k0)]]$concordance_dep_silhouette_color + 
              ggtitle("After filtering"))
      dev.off()
    }
    
    saveRDS(plots[[paste0("TPM_1_25p_", k0)]], 
            file = paste0(figdir, "/crossmethod_consistency_final", dtpext,
                          "_", k0, "_plots.rds"))
  }  
  concordances
}

summarize_runfailure <- function(figdir, datasets, exts, dtpext, cols,
                                 singledsfigdir, cobradir, concordancedir, 
                                 dschardir, origvsmockdir, distrdir, plotmethods, 
                                 dstypes, pch_ncells) {
  
  ## Initialize list to hold all plots
  plots <- list()
  
  ## Read all runstatus information
  runstatus <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/runfailure/", ds, e, 
                     "_runfailure_summary_data.rds"))$runstatus_sum %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods) %>%
        dplyr::mutate(dataset = paste0(dataset, "_", filt)) %>%
        dplyr::mutate(dataset = gsub("_$", "", dataset)) %>%
        dplyr::select(-filt)
    }))
  }))
  
  runstatus <- reshape2::dcast(runstatus, method ~ dataset, value.var = "failure_rate")
  rownames(runstatus) <- runstatus$method
  runstatus$method <- NULL

  ## --------------------------- Final summary plot ------------------------- ##
  pdf(paste0(figdir, "/runfailure_final", dtpext, ".pdf"), width = 10, height = 10)
  pheatmap::pheatmap(runstatus, cluster_rows = FALSE, cluster_cols = FALSE,
                     main = "Failure rates by data set", 
                     color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100))
  dev.off()

  runstatus
}
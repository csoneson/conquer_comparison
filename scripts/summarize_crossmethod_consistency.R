source("scripts/help_function_crossmethod_concordance.R")

summarize_crossmethod_consistency <- function(figdir, datasets, exts, dtpext, cols,
                                              singledsfigdir, cobradir, concordancedir, 
                                              dschardir, origvsmockdir) {
  pdf(paste0(figdir, "/summary_crossmethod_consistency", exts, dtpext, ".pdf"), 
      width = 14, height = 10)
  
  concordances <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      as.data.frame(readRDS(paste0(concordancedir, "/", ds, 
                                   e, "_concordances.rds"))$concordance_betweenmethods_pairwise)
    }))
  }))
  concordances$method1 <- gsub(paste(exts, collapse = "|"), "", concordances$method1)
  concordances$method2 <- gsub(paste(exts, collapse = "|"), "", concordances$method2)
  
  for (f in unique(concordances$filt)) {
    help_function_crossmethod_concordance(concordances %>% dplyr::filter(filt == f), 
                                          k0 = 100, titleext = f)
    help_function_crossmethod_concordance(concordances %>% dplyr::filter(filt == f), 
                                          k0 = 1000, titleext = f)
  }
  dev.off()
}
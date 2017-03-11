source("/home/Shared/data/seq/conquer/comparison/scripts/help_function_crossmethod_concordance.R")

summarize_crossmethod_consistency <- function(figdir, datasets, exts, dtpext, cols = cols) {
  pdf(paste0(figdir, "/summary_crossmethod_consistency", exts, dtpext, ".pdf"), width = 14, height = 10)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/consistency/", ds, exts, "_concordances.rds"))
  })
  concordances <- do.call(rbind, lapply(summary_data_list, 
                                        function(x) as.data.frame(x$concordance_betweenmethods_auc)))
  concordances$method1 <- gsub(exts, "", concordances$method1)
  concordances$method2 <- gsub(exts, "", concordances$method2)
  
  help_function_crossmethod_concordance(concordances)  
  dev.off()
}
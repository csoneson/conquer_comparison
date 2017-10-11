summarize_tsne <- function(figdir, datasets, exts, dtpext, cols,
                           singledsfigdir, cobradir, concordancedir, 
                           dschardir, origvsmockdir, distrdir, plotmethods, 
                           dstypes, pch_ncells) {
  
  X <- list()
  for (ds in datasets) {
    X[[ds]] <- readRDS(paste0(dschardir, "/", ds,
                              "_dataset_characteristics_plots.rds"))
  }
  pdf(paste0(figdir, "/tsne_final", dtpext, ".pdf"), width = 13.5, height = 4.25 * ceiling(length(X)/3))
  print(plot_grid(plotlist = lapply(
    X, function(x) {
      x$tsne + guides(color = guide_legend(nrow = 2))
    }), 
  ncol = 3, labels = LETTERS[seq_len(length(datasets))], align = "h"))
  dev.off()
  
  lapply(X, function(x) x$tsne)
}



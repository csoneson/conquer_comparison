#' Calculate "relative performance", using the results for the largest sample 
#' size as the truth. All performances (FDR, FPR, TPR) are calculated using each
#' of the methods' results as truth and represented in matrix form
#' 
plot_results_relativetruth_all <- function(cobra, relperf_alltruths, colvec, 
                                           exts, summary_data = list()) {

  ## Get unique set of sample sizes
  ttmp <- unique(as.numeric(get_nsamples(colnames(padj(cobra)))))
  ttmp <- as.character(sort(as.numeric(ttmp), decreasing = TRUE))
  for (m in ttmp) {
    for (k in unique(as.numeric(get_repl(colnames(padj(cobra)))))) {
      for (tb in names(relperf_alltruths)) {
        message(paste0(m, ".", k, ".", tb))
        tbt <- relperf_alltruths[[tb]]
        rownames(tbt) <- gsub(exts, "", rownames(tbt))
        colnames(tbt) <- gsub(exts, "", colnames(tbt))
        tbs <- get_nsamples(colnames(tbt)) == m & get_repl(colnames(tbt)) == k
        if (any(tbs)) {
          if (length(setdiff(unique(unlist(tbt[, tbs])), c(NaN, NA))) > 2) {
            pheatmap(tbt[, tbs], 
                     cluster_rows = FALSE, cluster_cols = FALSE, 
                     main = paste0(toupper(tb), ", all method truths (rows)"), display_numbers = TRUE,
                     annotation_col = data.frame(method = colnames(tbt[, tbs]), 
                                                 row.names = colnames(tbt[, tbs])), 
                     annotation_colors = list(method = structure(colvec, 
                                                                 names = paste0(gsub(exts, "", names(colvec)),
                                                                                ".", m, ".", k))),
                     annotation_legend = FALSE, annotation_names_col = FALSE)
          }
        }
      }
    }
  }
  return(invisible(summary_data))
}

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Calculate "relative performance", using the results for the largest sample 
#' size as the truth. All performances (FDR, FPR, TPR) are calculated using each
#' of the methods' results as truth and represented in matrix form
#' 
plot_results_relativetruth_all <- function(cobra, colvec, exts = exts, summary_data = list()) {
  fdr_all <- list()
  fpr_all <- list()
  tpr_all <- list()
  f1_all <- list()
  
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  for (m in gsub("\\.truth", "", grep("\\.truth", colnames(truth(cobratmp)), value = TRUE))) {
    cobrares <- calculate_performance(cobratmp, onlyshared = FALSE, 
                                      aspects = c("tpr", "fpr", "fdrtpr"), 
                                      binary_truth = paste0(m, ".truth"), 
                                      thrs = 0.05)
    fdr_all[[m]] <- data.frame(fdrtpr(cobrares)[, c("method", "FDR")], 
                               truth = paste0(m, " (", sum(truth(cobratmp)[, paste0(m, ".truth")], 
                                                           na.rm = TRUE), ")"),
                               stringsAsFactors = FALSE)
    fpr_all[[m]] <- data.frame(fpr(cobrares)[, c("method", "FPR")], 
                               truth = paste0(m, " (", sum(truth(cobratmp)[, paste0(m, ".truth")], 
                                                           na.rm = TRUE), ")"),
                               stringsAsFactors = FALSE)
    tpr_all[[m]] <- data.frame(tpr(cobrares)[, c("method", "TPR")], 
                               truth = paste0(m, " (", sum(truth(cobratmp)[, paste0(m, ".truth")], 
                                                           na.rm = TRUE), ")"),
                               stringsAsFactors = FALSE)
    tmp <- fdrtpr(cobrares) %>% dplyr::mutate(F1 = 2*TP/(2*TP + FN + FP))
    f1_all[[m]] <- data.frame(tmp[, c("method", "F1")],
                              truth = paste0(m, " (", sum(truth(cobratmp)[, paste0(m, ".truth")],
                                                          na.rm = TRUE), ")"),
                              stringsAsFactors = FALSE)
  }
  RES <- list(fdr = do.call(rbind, fdr_all) %>% dcast(truth ~ method, value.var = "FDR") %>% as.data.frame(),
              fpr = do.call(rbind, fpr_all) %>% dcast(truth ~ method, value.var = "FPR") %>% as.data.frame(),
              tpr = do.call(rbind, tpr_all) %>% dcast(truth ~ method, value.var = "TPR") %>% as.data.frame(),
              f1 = do.call(rbind, f1_all) %>% dcast(truth ~ method, value.var = "F1") %>% as.data.frame())
  RES <- lapply(RES, function(tb) {
    rownames(tb) <- tb$truth
    tb$truth <- NULL
    tb[order(rownames(tb)), ]
  })
  
  ttmp <- unique(as.numeric(get_nsamples(colnames(padj(cobratmp)))))
  ttmp <- as.character(sort(as.numeric(ttmp), decreasing = TRUE))
  for (m in ttmp) {
    for (k in unique(as.numeric(get_repl(colnames(padj(cobratmp)))))) {
      for (tb in names(RES)) {
        tbt <- RES[[tb]]
        tbs <- get_nsamples(colnames(tbt)) == m & get_repl(colnames(tbt)) == k
        if (any(tbs)) {
          if (length(setdiff(unique(unlist(tbt[, tbs])), c(NaN, NA))) > 2) {
            pheatmap(tbt[, tbs], 
                     cluster_rows = FALSE, cluster_cols = FALSE, 
                     main = paste0(toupper(tb), ", all method truths (rows)"), display_numbers = TRUE,
                     annotation_col = data.frame(method = colnames(tbt[, tbs]), 
                                                 row.names = colnames(tbt[, tbs])), 
                     annotation_colors = list(method = structure(colvec, names = paste0(names(colvec),
                                                                                        ".", m, ".", k))),
                     annotation_legend = FALSE, annotation_names_col = FALSE)
          }
        }
      }
    }
  }
  return(invisible(summary_data))
}

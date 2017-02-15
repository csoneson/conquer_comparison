source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Plot performance using each method's results with the largest sample size as
#' truth
#' 
plot_results_relativetruth <- function(cobra, colvec, summary_data = list()) {
  ## Generate a new cobradata object where the truth of each method is 
  ## considered to be the results obtained with the largest sample size.
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  tmp <- melt(as.matrix(padj(cobratmp))) %>% 
    tidyr::separate(Var2, into = c("method", "nsamples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(Var1 = paste0(method, ".", Var1)) 
  truth <- melt(as.matrix(truth(cobratmp)[, grep("\\.truth", colnames(truth(cobratmp)))])) %>%
    dplyr::mutate(gene = paste(gsub("\\.truth", "", Var2), Var1, sep = ".")) %>%
    dplyr::mutate(status = value) %>%
    dplyr::select(gene, status)
  rownames(truth) <- truth$gene
  truth$gene <- NULL
  ## Remove methods where the "truth" run failed
  tmp <- subset(tmp, method %in% unique(get_method(rownames(truth))))
  
  vals <- tmp %>% dplyr::select(Var1, Var2, value) %>% dcast(Var1 ~ Var2) %>% as.data.frame()
  rownames(vals) <- vals$Var1
  vals$Var1 <- NULL
  cobrarel <- COBRAData(padj = vals, truth = truth)
  
  cobrares <- calculate_performance(cobrarel, onlyshared = TRUE, 
                                    aspects = c("tpr", "fpr", "fdrtprcurve", "roc"), 
                                    binary_truth = "status", 
                                    thrs = 0.05)
  
  # ttmp <- unique(get_nsamples(basemethods(cobrares)))
  # ttmp <- as.character(sort(as.numeric(ttmp), decreasing = TRUE)[-1])
  # for (m in ttmp) {
  #   for (k in unique(get_repl(basemethods(cobrares)))) {
  #     c2 <- colvec
  #     names(c2) <- paste0(names(c2), ".", m, ".", k)
  #     km <- basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m & 
  #                                   get_repl(basemethods(cobrares)) == k]
  #     cobraplot <- prepare_data_for_plot(cobrares, keepmethods = km, colorscheme = c2[km])
  #     
  #     print(plot_fdrtprcurve(cobraplot, plottype = "curve") + ggtitle("Truth defined by each method"))
  #     print(plot_roc(cobraplot) + ggtitle("Truth defined by each method"))
  #   }
  # }
  
  ## Extract subset
  # for (m in unique(get_method(basemethods(cobrares)))) {
  #   plot_res_subset(cobrares, keepmethods = basemethods(cobrares)[get_method(basemethods(cobrares)) == m],
  #                   type = "method", colvec = colvec)
  # }
  
  for (m in unique(get_nsamples(basemethods(cobrares)))) {
    plot_res_subset(cobrares, keepmethods = basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m],
                    type = "number", colvec = colvec, nsamp = m)
  }
  return(invisible(summary_data))
}

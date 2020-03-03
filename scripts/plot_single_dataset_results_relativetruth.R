#' Help function to plot FPR and TPR for a subset of methods (defined either by
#' method or sample size)
plot_res_subset <- function(cobrares, keepmethods, type, colvec, nsamp = 1) {
  cobraplot <- prepare_data_for_plot(cobrares, keepmethods = keepmethods, 
                                     colorscheme = colvec)
  
  ## Modify method column so that all replicates with the same number of samples have the same name
  tpr(cobraplot) <- tpr(cobraplot) %>% 
    dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
  fpr(cobraplot) <- fpr(cobraplot) %>% 
    dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
  if (type == "method") {
    tmpvec <- colvec[seq_len(length(unique(tpr(cobraplot)$method)))]
    names(tmpvec) <- unique(tpr(cobraplot)$method)
  } else if (type == "number") {
    tmpvec <- colvec[sapply(strsplit(unique(tpr(cobraplot)$method), "\\."), .subset, 1)]
    names(tmpvec) <- paste0(names(tmpvec), ".", nsamp)
  }
  plotcolors(cobraplot) <- tmpvec
  print(plot_fpr(cobraplot, xaxisrange = c(0, min(1.1*max(fpr(cobrares)$FPR), 1))) + 
          ggtitle("Truth defined by each method"))
  print(plot_tpr(cobraplot) + ggtitle("Truth defined by each method"))
  
  print(ggplot(tpr(cobraplot), aes(x = DIFF, y = TPR, color = method)) + 
          geom_point(size = 5) + theme_bw() + xlab("Number of DE genes with maximal sample size") + 
          ylab("Relative TPR") + scale_color_manual(values = tmpvec) + 
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 13)))
  print(ggplot(fpr(cobraplot), aes(x = DIFF, y = FPR, color = method)) + 
          geom_point(size = 5) + theme_bw() + xlab("Number of DE genes with maximal sample size") + 
          ylab("Relative FPR") + scale_color_manual(values = tmpvec) + 
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 13)))
  print(ggplot(fpr(cobraplot), aes(x = DIFF, y = NBR, color = method)) + 
          geom_point(size = 5) + theme_bw() + xlab("Number of DE genes with maximal sample size") + 
          ylab("Number of significant genes") + scale_color_manual(values = tmpvec) + 
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 13)))
  
  return(list(fpr = fpr(cobraplot)[, c("thr", "basemethod", "FPR")], 
              tpr = tpr(cobraplot)[, c("thr", "basemethod", "TPR")]))
}


#' Plot performance using each method's results with the largest sample size as
#' truth
#' 
plot_results_relativetruth <- function(cobra, colvec, exts, summary_data = list()) {
  ## Generate a new cobradata object where the truth of each method is 
  ## considered to be the results obtained with the largest sample size.
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  ## Melt adjusted p-value matrix
  tmp <- melt(as.matrix(padj(cobratmp))) %>% 
    tidyr::separate(Var2, into = c("method", "nsamples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(Var1 = paste0(method, ".", Var1)) 
  ## Melt truth table
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
  ## Generate COBRAData object where genes are defined by method + gene
  cobrarel <- COBRAData(padj = vals, truth = truth)
  
  cobrares <- calculate_performance(cobrarel, onlyshared = TRUE, 
                                    aspects = c("tpr", "fpr"), 
                                    binary_truth = "status", 
                                    thrs = 0.05)
  
  for (m in unique(get_nsamples(basemethods(cobrares)))) {
    message(m)
    res <- plot_res_subset(cobrares, 
                           keepmethods = basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m],
                           type = "number", colvec = colvec, nsamp = m)
    
    summary_data$fpr_relative <- rbind(summary_data$fpr_relative, res$fpr)
    summary_data$tpr_relative <- rbind(summary_data$tpr_relative, res$tpr)
  }
  return(invisible(summary_data))
}

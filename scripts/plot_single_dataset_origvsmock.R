source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

plot_compare_orig_mock <- function(cobras, colvec, summary_data = list()) {
  ## Spearman correlations
  spearmans <- lapply(names(cobras), function(ncbr) {
    cbr <- cobras[[ncbr]]
    #pval(cbr)[is.na(pval(cbr))] <- 1
    #padj(cbr)[is.na(padj(cbr))] <- 1
    tmpmt <- pval(cbr)
    if (!(all(colnames(iCOBRA::score(cbr)) %in% colnames(tmpmt)))) {
      sdn <- setdiff(colnames(iCOBRA::score(cbr)), colnames(tmpmt))
      tmpmt <- merge(tmpmt, -iCOBRA::score(cbr)[, sdn, drop = FALSE], by = 0, all = TRUE)
      rownames(tmpmt) <- tmpmt$Row.names
      tmpmt$Row.names <- NULL
    }
    if (!(all(colnames(padj(cbr)) %in% colnames(tmpmt)))) {
      sdn <- setdiff(colnames(padj(cbr)), colnames(tmpmt))
      tmpmt <- merge(tmpmt, padj(cbr)[, sdn, drop = FALSE], by = 0, all = TRUE)
      rownames(tmpmt) <- tmpmt$Row.names
      tmpmt$Row.names <- NULL
    }
    tmpmt <- as.matrix(tmpmt)
    spcor <- cor(tmpmt, use = "pairwise.complete.obs", method = "spearman")
    spcor[is.na(spcor)] <- 0
    
    reshape2::melt(spcor) %>% 
      dplyr::rename(method1 = Var1) %>% 
      dplyr::rename(method2 = Var2) %>% 
      dplyr::filter(method1 != method2) %>% 
      tidyr::separate(method1, into = c("method1", "nbr_samples1", "replicate1"), sep = "\\.") %>%
      tidyr::separate(method2, into = c("method2", "nbr_samples2", "replicate2"), sep = "\\.") %>%
      dplyr::filter(method1 == method2 & nbr_samples1 == nbr_samples2) %>%
      dplyr::mutate(tp = ncbr) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_mock", "mock")) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_", "original"))
  })
  spm <- do.call(rbind, spearmans)
  
  summary_data$spearman <- spm
  
  for (nbrsamples in unique(intersect(subset(spm, tp == "original")$nbr_samples1,
                                      subset(spm, tp == "mock")$nbr_samples1))) {
    print(spm %>% dplyr::filter(nbr_samples1 == nbrsamples & nbr_samples2 == nbrsamples) %>%
            ggplot(aes(x = method1, y = value, color = method1, shape = tp)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + 
            ylab(paste0("Spearman correlation")) + 
            scale_color_manual(values = colvec, name = "method") + 
            scale_shape_discrete(name = "") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
            ggtitle(paste0(nbrsamples, " samples/group")))
  }
  
  ## Jaccard distances
  for (adjpthr in c(0.05, 0.1, 0.2)) {
    jaccards <- lapply(names(cobras), function(ncbr) {
      cbr <- cobras[[ncbr]]
      pval(cbr)[is.na(pval(cbr))] <- 1
      padj(cbr)[is.na(padj(cbr))] <- 1
      
      cobraperf <- calculate_performance(cbr, aspects = "overlap", 
                                         type_venn = "adjp", thr_venn = adjpthr)
      overlap(cobraperf) <- overlap(cobraperf)[, order(colnames(overlap(cobraperf)))]
      ol <- as.matrix(overlap(cobraperf))
      ol[is.na(ol)] <- 0
      
      jacc <- (t(ol) %*% ol)/(nrow(ol) - t(1 - ol) %*% (1 - ol))
      ## Set NaNs (caused by no significant genes) to 0
      w <- which(colSums(ol) == 0)
      for (i in w) {
        for (j in w) {
          if (is.na(jacc[i, j])) jacc[i, j] <- 0
        }
      }
      
      reshape2::melt(jacc) %>% 
        dplyr::rename(method1 = Var1) %>% 
        dplyr::rename(method2 = Var2) %>% 
        dplyr::filter(method1 != method2) %>% 
        tidyr::separate(method1, into = c("method1", "nbr_samples1", "replicate1"), sep = "\\.") %>%
        tidyr::separate(method2, into = c("method2", "nbr_samples2", "replicate2"), sep = "\\.") %>%
        dplyr::filter(method1 == method2 & nbr_samples1 == nbr_samples2) %>%
        dplyr::mutate(tp = ncbr) %>%
        dplyr::mutate(tp = replace(tp, tp == "tp_mock", "mock")) %>%
        dplyr::mutate(tp = replace(tp, tp == "tp_", "original"))
    })
    jaccm <- do.call(rbind, jaccards)
    
    summary_data[[paste0("jaccard_adjp", adjpthr)]] <- jaccm
    
    for (nbrsamples in unique(intersect(subset(jaccm, tp == "original")$nbr_samples1,
                                        subset(jaccm, tp == "mock")$nbr_samples1))) {
      print(jaccm %>% dplyr::filter(nbr_samples1 == nbrsamples & nbr_samples2 == nbrsamples) %>%
              ggplot(aes(x = method1, y = value, color = method1, shape = tp)) + 
              geom_point(size = 5) + theme_bw() + xlab("") + 
              ylab(paste0("Jaccard index (padj thr = ", adjpthr, ")")) + 
              scale_color_manual(values = colvec, name = "method") + 
              scale_shape_discrete(name = "") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              ggtitle(paste0(nbrsamples, " samples/group")))
    }
  }
  return(invisible(summary_data))
}

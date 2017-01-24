suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(lazyeval))
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")

#' Convenience functions to extract the first, second and third part of each of
#' a vector of strings (parts are separated by .)
#' 
get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

#' Make UpSet plot, making sure that the first and last columns have at least
#' one significant feature
#' 
plot_upset_with_reordering <- function(cobraplot, nintersects, ...) {
  ## Reorder so that the first and last columns have something significant
  m <- min(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m)) 
    overlap(cobraplot) <- overlap(cobraplot)[, c(m, setdiff(1:ncol(overlap(cobraplot)), m))]
  m <- max(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m))
    overlap(cobraplot) <- overlap(cobraplot)[, c(setdiff(1:ncol(overlap(cobraplot)), m), m)]
  tryCatch(plot_upset(cobraplot, order.by = "freq", decreasing = TRUE, nintersects = nintersects, ...),
           error = function(e) NULL)
}

#' Plot time used for each method
#' 
plot_timing <- function(timinglist, colvec, summary_data = list()) {
  timings <- sapply(timinglist, function(i) i["elapsed"])
  timings <- data.frame(method = names(timings), timing = timings) %>% 
    tidyr::separate(method, into = c("method", "nsamples", "repl", "elapsed"), sep = "\\.") %>%
    group_by(method, nsamples) %>% dplyr::summarise(timing = median(timing)) %>%
    dplyr::mutate(nsamples = as.numeric(as.character(nsamples)))
  summary_data$timing <- rbind(summary_data$timing, timings)
  print(timings %>%    
          ggplot(aes(x = nsamples, y = timing, group = method, color = method)) + 
          geom_line(size = 2.5) + scale_y_log10() + theme_bw() + 
          xlab("Number of replicates") + 
          scale_color_manual(values = colvec))
  return(invisible(summary_data))
}

#' Calculate "relative performance", using the results for the largest sample 
#' size as the truth. All performances (FDR, FPR, TPR) are calculated using each
#' of the methods' results as truth and represented in matrix form
#' 
plot_results_relativetruth_all <- function(cobra, colvec, summary_data = list()) {
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

#' Plot p-value and statistic for Kolmogorov-Smirnov test for uniformity of
#' p-values
#' 
plot_ks <- function(cobra, colvec, summary_data = list()) {
  ## For each method with p-values, calculate p-value from K-S test for uniformity
  pvs <- pval(cobra)
  ksp <- apply(pvs, 2, function(x) ks.test(x = x[!is.na(x)], y = punif, min = 0, max = 1)$p.value)
  kst <- apply(pvs, 2, function(x) ks.test(x = x[!is.na(x)], y = punif, min = 0, max = 1)$statistic)
  
  print(data.frame(method = names(kst), KSstat = kst, stringsAsFactors = FALSE) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = 
                          factor(nbr_samples, 
                                 levels = as.character(sort(as.numeric(as.character(nbr_samples)))))) %>% 
          ggplot(aes(x = method, y = KSstat, color = method, shape = nbr_samples)) + geom_point(size = 5) + 
          theme_bw() + xlab("") + ylab("Kolmogorov-Smirnov statistic, p-value distribution") + 
          scale_color_manual(values = colvec, name = "method") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  for (i in sort(unique(as.numeric(get_nsamples(names(kst)))))) {
    print(data.frame(method = names(kst), KSstat = kst, stringsAsFactors = FALSE) %>% 
            tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
            dplyr::filter(nbr_samples == i) %>% 
            ggplot(aes(x = method, y = KSstat, color = method)) + geom_point(size = 5) + 
            theme_bw() + xlab("") + ylab("Kolmogorov-Smirnov statistic, p-value distribution") + 
            scale_color_manual(values = colvec, name = "method") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(paste0(i, " samples per condition")))
  }
  return(invisible(summary_data))
}

#' Calculate and plot "true FPR" (fraction of genes with nominal p<0.05 in null
#' situation)
#' 
plot_truefpr <- function(cobra, colvec, summary_data = list()) {
  pvs <- pval(cobra)
  fpr <- apply(pvs, 2, function(x) length(which(x < 0.05))/length(x[!is.na(x)]))

  summary_data$fracpbelow0.05 <- 
    rbind(summary_data$fracpbelow0.05, 
          data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE))

  print(data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = 
                          factor(nbr_samples, 
                                 levels = as.character(sort(as.numeric(as.character(nbr_samples)))))) %>% 
          ggplot(aes(x = method, y = FPR, color = method, shape = nbr_samples)) + geom_point(size = 5) + 
          theme_bw() + xlab("") + ylab("Fraction of genes with nominal p < 0.05") + 
          geom_hline(yintercept = 0.05) + 
          scale_color_manual(values = colvec, name = "method") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  for (i in sort(unique(as.numeric(get_nsamples(names(fpr)))))) {
    print(data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE) %>% 
            tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
            dplyr::filter(nbr_samples == i) %>% 
            ggplot(aes(x = method, y = FPR, color = method)) + geom_point(size = 5) + 
            theme_bw() + xlab("") + ylab("Fraction of genes with nominal p < 0.05") + 
            geom_hline(yintercept = 0.05) + 
            scale_color_manual(values = colvec, name = "method") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(paste0(i, " samples per condition")))
  }
  return(invisible(summary_data))
}

#' Help function to plot FPR and TPR for a subset of methods (defined either by
#' method or sample size)
#' 
plot_res_subset <- function(cobrares, keepmethods, type, colvec, nsamp = 1) {
  cobraplot <- 
    prepare_data_for_plot(cobrares, keepmethods = keepmethods, colorscheme = colvec)
  
  ## Modify method column so that all replicates with the same number of samples have the same name
  tpr(cobraplot) <- tpr(cobraplot) %>% 
    dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
  fpr(cobraplot) <- fpr(cobraplot) %>% 
    dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
  if (type == "method") {
    tmpvec <- colvec[1:length(unique(tpr(cobraplot)$method))]
    names(tmpvec) <- unique(tpr(cobraplot)$method)
  } else if (type == "number") {
    tmpvec <- colvec[sapply(strsplit(unique(tpr(cobraplot)$method), "\\."), .subset, 1)]
    names(tmpvec) <- paste0(names(tmpvec), ".", nsamp)
  }
  plotcolors(cobraplot) <- tmpvec
  print(plot_fpr(cobraplot, xaxisrange = c(0, min(1.1*max(fpr(cobrares)$FPR), 1))) + 
          ggtitle("Truth defined by each method"))
  print(plot_tpr(cobraplot) + ggtitle("Truth defined by each method"))
}

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
  vals <- tmp %>% dplyr::select(Var1, Var2, value) %>% dcast(Var1 ~ Var2) %>% as.data.frame()
  rownames(vals) <- vals$Var1
  vals$Var1 <- NULL
  cobrarel <- COBRAData(padj = vals, truth = truth)
  
  cobrares <- calculate_performance(cobrarel, onlyshared = TRUE, 
                                    aspects = c("tpr", "fpr", "fdrtprcurve", "roc"), 
                                    binary_truth = "status", 
                                    thrs = 0.05)
  
  ttmp <- unique(get_nsamples(basemethods(cobrares)))
  ttmp <- as.character(sort(as.numeric(ttmp), decreasing = TRUE)[-1])
  for (m in ttmp) {
    for (k in unique(get_repl(basemethods(cobrares)))) {
      c2 <- colvec
      names(c2) <- paste0(names(c2), ".", m, ".", k)
      km <- basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m & 
                                    get_repl(basemethods(cobrares)) == k]
      cobraplot <- prepare_data_for_plot(cobrares, keepmethods = km, colorscheme = c2[km])
      
      print(plot_fdrtprcurve(cobraplot, plottype = "curve") + ggtitle("Truth defined by each method"))
      print(plot_roc(cobraplot) + ggtitle("Truth defined by each method"))
    }
  }
  
  ## Extract subset
  for (m in unique(get_method(basemethods(cobrares)))) {
    plot_res_subset(cobrares, keepmethods = basemethods(cobrares)[get_method(basemethods(cobrares)) == m],
                    type = "method", colvec = colvec)
  }
  
  for (m in unique(get_nsamples(basemethods(cobrares)))) {
    plot_res_subset(cobrares, keepmethods = basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m],
                    type = "number", colvec = colvec, nsamp = m)
  }
  return(invisible(summary_data))
}

#' Plot distribution of gene characteristics for significant/non-significant
#' genes found with each method
#' 
plot_results_characterization <- function(cobra, colvec, summary_data = list()) {
  sizes_nsamples <- gsub("avetpm.", "", grep("avetpm.", colnames(truth(cobra)), value = TRUE))
  for (szi in sizes_nsamples) {
    message(szi)
    pvals <- pval(cobra)[, grep(paste0("\\.", szi, "$"), colnames(pval(cobra)))]
    pvals[is.na(pvals)] <- 1
    padjs <- padj(cobra)[, grep(paste0("\\.", szi, "$"), colnames(padj(cobra)))]
    padjs[is.na(padjs)] <- 1
    tth <- truth(cobra)[, grep(paste0("\\.", szi, "$"), colnames(truth(cobra)))]
    colnames(tth) <- gsub(paste0("\\.", szi), "", colnames(tth))

    ## Number of methods calling each variable significant
    tmp1 <- rowSums(padjs[match(rownames(tth), rownames(padjs)), ] <= 0.05)
    tth$nbr_methods <- tmp1[match(rownames(tth), names(tmp1))]
    
    df2 <- merge(melt(as.matrix(padjs)), tth, by.x = "Var1", by.y = 0, all = TRUE)
    df2 <- subset(df2, rowSums(is.na(df2)) == 0)
    df2$sign <- df2$value <= 0.05
    df2$Var2 <- factor(df2$Var2, levels = sort(levels(df2$Var2)))
    col_reorder <- colvec
    names(col_reorder) <- paste0(names(col_reorder), ".", szi)
    for (y in c("avetpm", "avecount", "vartpm")) {
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~log2(x), x = as.name(y))) %>%
                group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch((median(newy[sign]) - median(newy[!sign]))/median(newy),
                                        error = function(e) NA)) %>%
                dplyr::mutate(charac = paste0("log2_", y)))

      print(ggplot(df2, aes_string(x = "Var2", y = paste0("log2(", y, ")"), 
                                   fill = "Var2", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = log2(max(df2[, y])),
                                                            label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              ylab(paste0(expression(log[2]), "(", ifelse(y == "avetpm", "average TPM)", 
                                                          ifelse(y == "avecount", "average count)", 
                                                                 "variance(TPM))")))))
      print(ggplot(df2, aes_string(x = paste0("log2(", y, ")"), fill = "Var2", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var2, scales = "free_y") + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(paste0(expression(log[2]), "(", ifelse(y == "avetpm", "average TPM)", 
                                                          ifelse(y == "avecount", "average count)", 
                                                                 "variance(TPM))")))))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = paste0("log2(", y, ")"), 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + ylab(paste0(expression(log[2]), "(",
                                                        ifelse(y == "avetpm", "average TPM)", 
                                                               ifelse(y == "avecount", "average count)",
                                                                      "variance(TPM))")))) + 
              xlab("Number of methods calling gene significant"))
    }
    for (y in c("fraczero", "fraczerodiff")) {
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~x, x = as.name(y))) %>%
                group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch((median(newy[sign]) - median(newy[!sign]))/median(newy),
                                        error = function(e) NA)) %>%
                dplyr::mutate(charac = y))
      
      print(ggplot(df2, aes_string(x = "Var2", y = y, fill = "Var2", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = max(df2[, y]), label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              ylab(ifelse(y == "fraczero", "Zero fraction", "Difference in zero fraction")))
      print(ggplot(df2, aes_string(x = y, fill = "Var2", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var2, scales = "free_y") + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(ifelse(y == "fraczero", "Zero fraction", "Difference in zero fraction")))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = y, 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + ylab(ifelse(y == "fraczero", "Zero fraction", 
                                                        "Difference in zero fraction")) + 
              xlab("Number of methods calling gene significant"))
    }
  }
  return(invisible(summary_data))
}

compare_orig_mock <- function(cobras, colvec, summary_data = list()) {
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
  
  for (nbrsamples in unique(intersect(subset(spm, tp == "original")$nbr_samples1,
                                      subset(spm, tp == "mock")$nbr_samples1))) {
    print(spm %>% dplyr::filter(nbr_samples1 == nbrsamples & nbr_samples2 == nbrsamples) %>%
            ggplot(aes(x = method1, y = value, color = method1, alpha = tp)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + 
            ylab(paste0("Spearman correlation")) + 
            scale_color_manual(values = colvec, name = "method") + 
            scale_alpha_discrete(name = "") + 
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
    
    for (nbrsamples in unique(intersect(subset(jaccm, tp == "original")$nbr_samples1,
                                        subset(jaccm, tp == "mock")$nbr_samples1))) {
      print(jaccm %>% dplyr::filter(nbr_samples1 == nbrsamples & nbr_samples2 == nbrsamples) %>%
              ggplot(aes(x = method1, y = value, color = method1, alpha = tp)) + 
              geom_point(size = 5) + theme_bw() + xlab("") + 
              ylab(paste0("Jaccard index (padj thr = ", adjpthr, ")")) + 
              scale_color_manual(values = colvec, name = "method") + 
              scale_alpha_discrete(name = "") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              ggtitle(paste0(nbrsamples, " samples/group")))
    }
  }
  return(invisible(summary_data))
}

plot_results <- function(cobra, colvec, summary_data = list()) {
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  cobraperf <- calculate_performance(cobratmp, aspects = "overlap", 
                                     type_venn = "adjp", thr_venn = 0.05)
  overlap(cobraperf) <- overlap(cobraperf)[, order(colnames(overlap(cobraperf)))]
  ol <- as.matrix(overlap(cobraperf))
  ol[is.na(ol)] <- 0
  
  nosignif <- colnames(ol)[which(colSums(ol) == 0)]
  
  ## --------------------- Number of detections ----------------------------- ##
  print(reshape2::melt(ol, varnames = c("gene", "method")) %>% 
          group_by(method) %>% 
          dplyr::summarise(nsignif = sum(value)) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = factor(nbr_samples, levels = as.character(unique(sort(as.numeric(as.character(nbr_samples))))))) %>%
          ggplot(aes(x = method, y = nsignif, color = method, shape = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + xlab("") + ylab("Number of detections") + 
          scale_color_manual(values = colvec) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  ## Divide by total number of genes
  print(reshape2::melt(sweep(ol, 2, colSums(!is.na(ol)), "/"), 
                       varnames = c("gene", "method")) %>% 
          group_by(method) %>% 
          dplyr::summarise(nsignif = sum(value)) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = factor(nbr_samples, levels = as.character(unique(sort(as.numeric(as.character(nbr_samples))))))) %>%
          ggplot(aes(x = method, y = nsignif, color = method, shape = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + xlab("") + ylab("Fraction of genes called significant") + 
          scale_color_manual(values = colvec) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  ## ----------------------- Jaccard distance ------------------------------- ##
  jacc <- (t(ol) %*% ol)/(nrow(ol) - t(1 - ol) %*% (1 - ol))
  ## Set NaNs (caused by no significant genes) to 0
  w <- which(colSums(ol) == 0)
  for (i in w) {
    for (j in w) {
      if (is.na(jacc[i, j])) jacc[i, j] <- 0
    }
  }
  
  mdsjacc <- cmdscale(1 - jacc, k = 2, eig = TRUE)
  mdsjaccx <- data.frame(mdsjacc$points)
  colnames(mdsjaccx) <- c("MDS_Jaccard_1", "MDS_Jaccard_2")
  mdsjaccx$mth <- rownames(mdsjaccx)
  
  ## Only methods actually detecting something
  jacc2 <- jacc[-which(rownames(jacc) %in% nosignif), -which(colnames(jacc) %in% nosignif)]
  mdsjacc2 <- cmdscale(1 - jacc2, k = 2, eig = TRUE)
  mdsjaccx2 <- data.frame(mdsjacc2$points)
  colnames(mdsjaccx2) <- c("MDS_Jaccard_1", "MDS_Jaccard_2")
  mdsjaccx2$mth <- rownames(mdsjaccx2)
  
  ## ---------------------- Spearman correlations --------------------------- ##
  tmpmt <- pval(cobratmp)
  if (!(all(colnames(iCOBRA::score(cobratmp)) %in% colnames(tmpmt)))) {
    sdn <- setdiff(colnames(iCOBRA::score(cobratmp)), colnames(tmpmt))
    tmpmt <- merge(tmpmt, -iCOBRA::score(cobratmp)[, sdn, drop = FALSE], by = 0, all = TRUE)
    rownames(tmpmt) <- tmpmt$Row.names
    tmpmt$Row.names <- NULL
  }
  if (!(all(colnames(padj(cobratmp)) %in% colnames(tmpmt)))) {
    sdn <- setdiff(colnames(padj(cobratmp)), colnames(tmpmt))
    tmpmt <- merge(tmpmt, padj(cobratmp)[, sdn, drop = FALSE], by = 0, all = TRUE)
    rownames(tmpmt) <- tmpmt$Row.names
    tmpmt$Row.names <- NULL
  }
  # for (i in 1:ncol(tmpmt)) {
  #   tmpmt[which(is.na(tmpmt[, i])), i] <- max(tmpmt[which(!is.na(tmpmt[, i])), i])
  # }
  tmpmt <- as.matrix(tmpmt)
  spdist <- 1 - cor(tmpmt, use = "pairwise.complete.obs", method = "spearman")
  spdist[is.na(spdist)] <- 2
  mdssp <- cmdscale(spdist, k = 2, eig = TRUE)
  mdsspx <- data.frame(mdssp$points)
  colnames(mdsspx) <- c("MDS_Spearman_1", "MDS_Spearman_2")
  mdsspx$mth <- rownames(mdsspx)
  
  ## ------------------------ Get summary info ------------------------------ ##
  tmpinfo <- data.frame(mth = mdsjaccx$mth) %>% 
    tidyr::separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.")
  
  all_methods <- unique(tmpinfo$method)
  all_sizes <- as.numeric(unique(tmpinfo$nbr_samples))
  all_replicates <- as.numeric(unique(tmpinfo$replicate))
  max_nreps <- max(as.numeric(tmpinfo$replicate))
  
  ## Only methods detecting anything
  tmpinfo2 <- data.frame(mth = mdsjaccx2$mth) %>% 
    tidyr::separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.")
  
  all_methods2 <- unique(tmpinfo2$method)
  all_sizes2 <- as.numeric(unique(tmpinfo2$nbr_samples))
  all_replicates2 <- as.numeric(unique(tmpinfo2$replicate))
  max_nreps2 <- max(as.numeric(tmpinfo2$replicate))
  
  ## ------------------ Number of detections shared between methods --------- ##
  nshared <- lapply(structure(all_sizes, names = all_sizes), function(i) {
    lapply(structure(all_replicates, names = all_replicates), function(rpl) {
      tmpol <- ol[, get_nsamples(colnames(ol)) == i & get_repl(colnames(ol)) == rpl]
      if (ncol(tmpol) > 0) {
        rsms <- rowSums(tmpol)
        rsms <- rsms[rsms != 0]
        rsms <- table(rsms)
        rev(cumsum(rev(rsms)))
      } else {
        NULL
      }
    })
  })
  tmp1 <- unlist(nshared, use.names = TRUE)
  df <- data.frame(group = names(tmp1), value = tmp1) %>%
    tidyr::separate(group, into = c("nsamples", "replicate", "nmethods"), sep = "\\.")
  for (i in all_sizes) {
    for (j in all_replicates) {
      s <- subset(df, nsamples == i & replicate == j)
      if (nrow(s) > 0) {
        msn <- setdiff(1:length(all_methods), s$nmethods)
        if (length(msn) > 0) {
          df <- rbind(df, data.frame(nsamples = i, replicate = j, nmethods = msn, value = 0))
        }
      }
    }
  }
  df$nmethods <- as.numeric(df$nmethods)
  df$nsamples <- factor(df$nsamples, 
                        levels = as.character(unique(sort(as.numeric(as.character(df$nsamples))))))
  ## Actual number
  print(
    df %>% 
      ggplot(aes(x = nmethods, y = value, 
                 group = interaction(nsamples, replicate), color = nsamples)) + 
      geom_line(size = 2) + theme_bw() + 
      scale_color_manual(values = structure(colvec[c(1,4,8,7,2,6,5,10,3,9,11,12)], names = NULL), 
                         name = "Number of\ncells\nper group") + 
      xlab("Number of methods") + ylab("Number of genes shared by each number of methods") + 
      scale_x_continuous(breaks = 1:length(all_methods), labels = 1:length(all_methods)) + 
      theme(legend.justification = c(1, 1), legend.position = c(1, 1)))
  
  ## Fraction
  print(
    df %>% group_by(nsamples, replicate) %>% 
      dplyr::mutate(value = value/max(value)) %>% 
      ggplot(aes(x = nmethods, y = value, 
                 group = interaction(nsamples, replicate), color = nsamples)) + 
      geom_line(size = 2) + theme_bw() + 
      scale_color_manual(values = structure(colvec[c(1,4,8,7,2,6,5,10,3,9,11,12)], names = NULL), 
                         name = "Number of\ncells\nper group") + 
      xlab("Number of methods") + 
      ylab("Fraction of all detected genes shared by each number of methods") + 
      scale_x_continuous(breaks = 1:length(all_methods), labels = 1:length(all_methods)) + 
      theme(legend.justification = c(1, 1), legend.position = c(1, 1)))
  
  ## ------------------ Jaccard coefficient distributions ------------------- ##
  jaccm <- reshape2::melt(jacc) %>% 
    dplyr::rename(method1 = Var1) %>% 
    dplyr::rename(method2 = Var2) %>% 
    dplyr::filter(method1 != method2) %>% 
    tidyr::separate(method1, into = c("method1", "nbr_samples1", "replicate1"), sep = "\\.") %>%
    tidyr::separate(method2, into = c("method2", "nbr_samples2", "replicate2"), sep = "\\.") %>%
    dplyr::filter(method1 == method2)
  
  print(jaccm %>%
          dplyr::select(method1, value) %>%
          ggplot(aes(x = method1, y = value, color = method1)) + 
          geom_point(size = 5) + theme_bw() + xlab("") + ylab("Jaccard index") + 
          scale_color_manual(values = colvec, name = "method") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))  
  
  for (sz in all_sizes) {
    print(jaccm %>%
            dplyr::filter(nbr_samples1 == sz) %>%
            dplyr::filter(nbr_samples2 == sz) %>%
            ggplot(aes(x = method1, y = value, color = method1)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + ylab("Jaccard index") + 
            scale_color_manual(values = colvec, name = "method") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
            ggtitle(paste0(sz, " samples/group")))
  }
  
  ## ----------------- Spearman correlation distributions ------------------- ##
  spdistm <- reshape2::melt(1 - spdist) %>% 
    dplyr::rename(method1 = Var1) %>% 
    dplyr::rename(method2 = Var2) %>% 
    dplyr::filter(method1 != method2) %>% 
    tidyr::separate(method1, into = c("method1", "nbr_samples1", "replicate1"), sep = "\\.") %>%
    tidyr::separate(method2, into = c("method2", "nbr_samples2", "replicate2"), sep = "\\.") %>%
    dplyr::filter(method1 == method2)
  
  print(spdistm %>%
          dplyr::select(method1, value) %>%
          ggplot(aes(x = method1, y = value, color = method1)) + 
          geom_point(size = 5) + theme_bw() + xlab("") + ylab("Spearman correlation") + 
          scale_color_manual(values = colvec, name = "method") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))  
  
  for (sz in all_sizes) {
    print(spdistm %>%
            dplyr::filter(nbr_samples1 == sz) %>%
            dplyr::filter(nbr_samples2 == sz) %>%
            ggplot(aes(x = method1, y = value, color = method1)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + ylab("Spearman correlation") + 
            scale_color_manual(values = colvec, name = "method") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
            ggtitle(paste0(sz, " samples/group")))
  }
  
  ## ------------------------- Plot global MDS ------------------------------ ##
  print(mdsjaccx %>% 
          separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
          ggplot(aes(x = MDS_Jaccard_1, y = MDS_Jaccard_2, col = method, pch = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + scale_color_manual(values = colvec))
  print(data.frame(dimension = 1:length(mdsjacc$eig),
                   eigenvalue = mdsjacc$eig) %>% 
          ggplot(aes(x = dimension, y = eigenvalue)) + geom_line() + 
          geom_point(size = 5) + ggtitle("MDS_Jaccard"))
  print(mdsjaccx2 %>% 
          separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
          ggplot(aes(x = MDS_Jaccard_1, y = MDS_Jaccard_2, col = method, pch = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + scale_color_manual(values = colvec))
  print(data.frame(dimension = 1:length(mdsjacc2$eig),
                   eigenvalue = mdsjacc2$eig) %>% 
          ggplot(aes(x = dimension, y = eigenvalue)) + geom_line() + 
          geom_point(size = 5) + ggtitle("MDS_Jaccard"))
  
  print(mdsspx %>% 
          separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
          ggplot(aes(x = MDS_Spearman_1, y = MDS_Spearman_2, col = method, pch = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + scale_color_manual(values = colvec))
  print(data.frame(dimension = 1:length(mdssp$eig),
                   eigenvalue = mdssp$eig) %>% 
          ggplot(aes(x = dimension, y = eigenvalue)) + geom_line() + 
          geom_point(size = 5) + ggtitle("MDS_Spearman"))
  if (length(all_methods) > 1) {
    ## MDS
    for (sz in all_sizes) {
      for (j in all_replicates) {
        ## A specific dataset (number of samples/replicate combination)
        keepmethods <- intersect(paste0(all_methods, ".", sz, ".", j),
                                 rownames(spdist))
        if (length(keepmethods) > 0) {
          mdssub <- cmdscale(spdist[keepmethods, keepmethods], k = 2, eig = TRUE)
          mdssubx <- data.frame(mdssub$points)
          colnames(mdssubx) <- c("MDS_Spearman_1", "MDS_Spearman_2")
          mdssubx$mth <- rownames(mdssubx)
          print(mdssubx %>% 
                  separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
                  ggplot(aes(x = MDS_Spearman_1, y = MDS_Spearman_2, 
                             col = method, pch = nbr_samples)) + 
                  geom_point(size = 10) + theme_bw() + 
                  scale_color_manual(values = colvec, guide = FALSE) + 
                  geom_label_repel(aes(label = method), size = 7) + 
                  theme(legend.position = "bottom"))
        }
      }
      
      ## All replicates, given number of samples
      keepmethods <- intersect(unlist(lapply(paste0(all_methods, ".", sz), 
                                             function(m) paste0(m, ".", all_replicates))), 
                               rownames(spdist))
      if (length(keepmethods) > 0) {
        mdssub <- cmdscale(spdist[keepmethods, keepmethods], k = 2, eig = TRUE)
        mdssubx <- data.frame(mdssub$points)
        colnames(mdssubx) <- c("MDS_Spearman_1", "MDS_Spearman_2")
        mdssubx$mth <- rownames(mdssubx)
        print(mdssubx %>% 
                separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
                ggplot(aes(x = MDS_Spearman_1, y = MDS_Spearman_2, 
                           col = method, pch = nbr_samples)) + 
                geom_point(size = 5) + theme_bw() + 
                scale_color_manual(values = colvec, guide = FALSE) + 
                geom_label_repel(aes(label = method)) + 
                theme(legend.position = "bottom"))
      }
    }
    
    ## Hierarchical clustering and heatmap
    for (sz in all_sizes) {
      for (j in all_replicates) {
        keepmethods <- intersect(paste0(all_methods, ".", sz, ".", j), 
                                 rownames(spdist))
        if (length(keepmethods) > 0) {
          hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
          plot(hcl)
          hclspd <- 1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order]
          pheatmap(hclspd, 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   scale = "none", main = "Spearman correlation", display_numbers = TRUE,
                   annotation_col = data.frame(method = colnames(hclspd), 
                                               row.names = colnames(hclspd)), 
                   annotation_colors = list(method = structure(colvec, names = paste0(names(colvec),
                                                                                      ".", sz, ".", j))),
                   annotation_legend = FALSE, annotation_names_col = FALSE,
                   annotation_row = data.frame(method = rownames(hclspd), 
                                               row.names = rownames(hclspd)),
                   annotation_names_row = FALSE)
        }
      }
      ## All replicates
      keepmethods <- intersect(unlist(lapply(paste0(all_methods, ".", sz), 
                                             function(m) paste0(m, ".", all_replicates))), 
                               rownames(spdist))
      if (length(keepmethods) > 0) {
        hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
        plot(hcl)
        hclspd <- 1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order]
        pheatmap(hclspd, 
                 cluster_rows = FALSE, cluster_cols = FALSE, 
                 scale = "none", main = "Spearman correlation", display_numbers = TRUE,
                 annotation_col = data.frame(method = colnames(hclspd), 
                                             row.names = colnames(hclspd)), 
                 annotation_colors = 
                   list(method = structure(rep(colvec, each = length(all_replicates)), 
                                           names = unlist(lapply(paste0(names(colvec), ".", sz), 
                                                                 function(m) paste0(m, ".", 
                                                                                    all_replicates))))),
                 annotation_legend = FALSE, annotation_names_col = FALSE,
                 annotation_row = data.frame(method = rownames(hclspd), 
                                             row.names = rownames(hclspd)),
                 annotation_names_row = FALSE)
      }
    }
    
    ## UpSet plots
    for (sz in all_sizes) {
      for (j in all_replicates) {
        c2 <- colvec
        names(c2) <- paste0(names(c2), ".", sz, ".", j)
        km <- paste0(all_methods, ".", sz, ".", j)
        cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
                                     colorscheme = c2[km], incloverall = FALSE)
        if (ncol(overlap(cpl)) > 0) {
          plot_upset_with_reordering(cpl, nintersects = 25,
                                     set.metadata = list(data = data.frame(sets = km,
                                                                           mth = km,
                                                                           row.names = km),
                                                         plots = list(list(type = "matrix_rows",
                                                                           column = "mth",
                                                                           colors = c2[km], alpha = 0.25))))
        }
      }
    }
    
    ## Include attribute plot    
    for (sz in all_sizes) {
      for (j in all_replicates) {
        c2 <- colvec
        names(c2) <- paste0(names(c2), ".", sz, ".", j)
        km <- paste0(all_methods, ".", sz, ".", j)
        cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
                                     colorscheme = c2[km], incloverall = FALSE)
        if (ncol(overlap(cpl)) > 0) {
          overlap_table <- overlap(cpl)
          m <- min(which(colSums(overlap_table) > 0))
          if (is.finite(m)) 
            overlap_table <- overlap_table[, c(m, setdiff(1:ncol(overlap_table), m))]
          m <- max(which(colSums(overlap_table) > 0))
          if (is.finite(m))
            overlap_table <- overlap_table[, c(setdiff(1:ncol(overlap_table), m), m)]
          plotorder <- colnames(overlap_table)[order(colSums(overlap_table), 
                                                     seq(1:ncol(overlap_table)),
                                                     decreasing = "true")]
          nsets <- ncol(overlap_table)
          nintersects <- 25
          sets.bar.color <- plotcolors(cpl)[plotorder]
          overlap_table$fraczero <- truth(cobratmp)[match(rownames(overlap_table), 
                                                          rownames(truth(cobratmp))), 
                                                 paste0("fraczero.", sz, ".", j)]
          upset(overlap_table, nsets = nsets, nintersects = nintersects, 
                sets.bar.color = sets.bar.color, 
                order.by = "freq", decreasing = TRUE, boxplot.summary = "fraczero",
                set.metadata = list(data = data.frame(sets = km, 
                                                      mth = km,
                                                      row.names = km),
                                    plots = list(list(type = "matrix_rows", 
                                                      column = "mth", 
                                                      colors = c2[km], alpha = 0.25))))
        }
      }
    }
    
    for (sz in all_sizes) {
      for (j in all_replicates) {
        c2 <- colvec
        names(c2) <- paste0(names(c2), ".", sz, ".", j)
        km <- paste0(all_methods, ".", sz, ".", j)
        cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
                                     colorscheme = c2[km], incloverall = FALSE)
        if (ncol(overlap(cpl)) > 0) {
          for (k in all_methods) {
            tmp_cp <- as(cpl, "COBRAPerformance")
            overlap(tmp_cp) <- overlap(tmp_cp)[overlap(tmp_cp)[, paste0(k, ".", sz, ".", j)] == 0 & 
                                                 rowSums(overlap(tmp_cp)) > 0, ]
            km <- paste0(setdiff(all_methods, k), ".", sz, ".", j)
            cpl1 <- prepare_data_for_plot(tmp_cp, keepmethods = km, 
                                          colorscheme = c2[km], incloverall = FALSE)
            plot_upset_with_reordering(cpl1, nintersects = 25, 
                                       mainbar.y.label = paste0("Intersection Size (only genes ", 
                                                                "not called DE by ", k, ")"), 
                                       set.metadata = list(data = data.frame(sets = km, 
                                                                             mth = km,
                                                                             row.names = km),
                                                           plots = list(list(type = "matrix_rows", 
                                                                             column = "mth", 
                                                                             colors = c2[km], alpha = 0.25))))
          }
        }
      }
    }
  }
  if (length(all_sizes) > 1) {
    for (mth in all_methods) {
      c2 <- colvec[mth]
      km <- intersect(unlist(lapply(paste0(mth, ".", all_sizes), 
                                    function(m) paste0(m, ".", all_replicates))), 
                      basemethods(cobraperf))
      cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
                                   colorscheme = structure(rep(c2, length(km)), names = km), 
                                   incloverall = FALSE)
      if (ncol(overlap(cpl)) > 0)
        plot_upset_with_reordering(cpl, nintersects = 25, 
                                   set.metadata = list(data = data.frame(sets = km, 
                                                                         mth = km,
                                                                         row.names = km),
                                                       plots = list(list(type = "matrix_rows", 
                                                                         column = "mth", 
                                                                         colors = c2[km], alpha = 0.25))))
    }
  }
  if (max_nreps > 1) {
    for (mth in all_methods) {
      for (sz in all_sizes) {
        c2 <- colvec[mth]
        tmpmth <- intersect(paste0(mth, ".", sz, ".", all_replicates), basemethods(cobraperf))
        if (length(tmpmth) > 1) {
          cpl <- prepare_data_for_plot(cobraperf, keepmethods = tmpmth, 
                                       colorscheme = structure(rep(c2, length(tmpmth)), names = tmpmth), 
                                       incloverall = FALSE)
          if (ncol(overlap(cpl)) > 0)
            plot_upset_with_reordering(cpl, nintersects = 25, 
                                       set.metadata = list(data = data.frame(sets = tmpmth, 
                                                                             mth = tmpmth,
                                                                             row.names = tmpmth),
                                                           plots = list(list(type = "matrix_rows", 
                                                                             column = "mth", 
                                                                             colors = c2[tmpmth], 
                                                                             alpha = 0.25))))
        }
      }
    }
  }
  return(invisible(summary_data))
}


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
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

plot_upset_with_reordering <- function(cobraplot, nintersects) {
  m <- min(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m)) 
    overlap(cobraplot) <- overlap(cobraplot)[, c(m, setdiff(1:ncol(overlap(cobraplot)), m))]
  m <- max(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m))
    overlap(cobraplot) <- overlap(cobraplot)[, c(setdiff(1:ncol(overlap(cobraplot)), m), m)]
  tryCatch(plot_upset(cobraplot, order.by = "freq", decreasing = TRUE, nintersects = nintersects),
           error = function(e) NULL)
}

plot_timing <- function(timinglist, colvec) {
  timings <- sapply(timinglist, function(i) i["elapsed"])
  timings <- data.frame(method = names(timings), timing = timings) %>% 
    tidyr::separate(method, into = c("method", "nsamples", "repl", "elapsed")) %>%
    group_by(method, nsamples) %>% dplyr::summarise(timing = median(timing)) %>%
    dplyr::mutate(nsamples = as.numeric(as.character(nsamples)))
  print(timings %>%    
          ggplot(aes(x = nsamples, y = timing, group = method, color = method)) + 
          geom_line(size = 2.5) + geom_point() + scale_y_log10() + theme_bw() + 
          scale_color_manual(values = colvec))
}

plot_results_relativetruth <- function(cobra, colvec) {
  ## Generate a new cobradata object where the truth of each method is 
  ## considered to be the results obtained with the largest sample size.
  maxn <- max(as.numeric(sapply(strsplit(colnames(padj(cobra)), "\\."), .subset, 2)))
  tmp <- melt(as.matrix(padj(cobra))) %>% 
    tidyr::separate(Var2, into = c("method", "nsamples", "repl"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(Var1 = paste0(method, ".", Var1)) 
  truth <- tmp %>%
    dplyr::filter(nsamples == maxn) %>% 
    dplyr::mutate(status = as.numeric(value < 0.05)) %>%
    dplyr::rename(gene = Var1) %>% 
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
      cobraplot <- 
        prepare_data_for_plot(cobrares, 
                              keepmethods = basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m & 
                                                                    get_repl(basemethods(cobrares)) == k], 
                              colorscheme = "Set3")
      
      ## Modify method column so that all replicates with the same number of samples have the same name
      print(plot_fdrtprcurve(cobraplot, plottype = "curve"))
      print(plot_roc(cobraplot))
    }
  }
  
  ## Extract subset
  for (m in unique(get_method(basemethods(cobrares)))) {
    cobraplot <- 
      prepare_data_for_plot(cobrares, 
                            keepmethods = basemethods(cobrares)[get_method(basemethods(cobrares)) == m], 
                            colorscheme = "Set3")
    
    ## Modify method column so that all replicates with the same number of samples have the same name
    tpr(cobraplot) <- tpr(cobraplot) %>% 
      dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
    fpr(cobraplot) <- fpr(cobraplot) %>% 
      dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
    tmpvec <- colvec[1:length(unique(tpr(cobraplot)$method))]
    names(tmpvec) <- unique(tpr(cobraplot)$method)
    plotcolors(cobraplot) <- tmpvec
    print(plot_fpr(cobraplot, xaxisrange = c(0, min(1.1*max(fpr(cobrares)$FPR), 1))))
    print(plot_tpr(cobraplot))
  }
  
  for (m in unique(get_nsamples(basemethods(cobrares)))) {
    cobraplot <- 
      prepare_data_for_plot(cobrares, 
                            keepmethods = basemethods(cobrares)[get_nsamples(basemethods(cobrares)) == m], 
                            colorscheme = "Set3")
    
    ## Modify method column so that all replicates with the same number of samples have the same name
    tpr(cobraplot) <- tpr(cobraplot) %>% 
      dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
    fpr(cobraplot) <- fpr(cobraplot) %>% 
      dplyr::mutate(method = paste0(get_method(method), ".", get_nsamples(method)))
    tmpvec <- colvec[1:length(unique(tpr(cobraplot)$method))]
    names(tmpvec) <- unique(tpr(cobraplot)$method)
    plotcolors(cobraplot) <- tmpvec
    print(plot_fpr(cobraplot, xaxisrange = c(0, min(1.1*max(fpr(cobrares)$FPR), 1))))
    print(plot_tpr(cobraplot))
  }
}

plot_results_characterization <- function(cobra, config, colvec) {
  ## TODO: Add some measure of outliers/non-homogeneity/multimodality?
  config <- fromJSON(file = config)
  mae <- readRDS(config$mae)
  groupid <- config$groupid
  mae <- clean_mae(mae = mae, groupid = groupid)
  
  subsets <- readRDS(config$subfile)
  keep_samples <- subsets$keep_samples
  imposed_condition <- subsets$out_condition
  
  sizes <- names(keep_samples)
  for (sz in sizes) {
    for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
      message(sz, ".", i)
      L <- subset_mae(mae, keep_samples, sz, i, imposed_condition)
      pvals <- pval(cobra)[, grep(paste0("\\.", sz, "\\.", i, "$"), colnames(pval(cobra)))]
      pvals[is.na(pvals)] <- 1
      padjs <- padj(cobra)[, grep(paste0("\\.", sz, "\\.", i, "$"), colnames(padj(cobra)))]
      padjs[is.na(padjs)] <- 1
      avecount <- data.frame(avecount = apply(L$count, 1, mean))
      avetpm <- data.frame(avetpm = apply(L$tpm, 1, mean))
      fraczero <- data.frame(fraczero = apply(L$count, 1, function(x) mean(x == 0)),
                             fraczero1 = apply(L$count[, L$condt == levels(factor(L$condt))[1]], 
                                               1, function(x) mean(x == 0)),
                             fraczero2 = apply(L$count[, L$condt == levels(factor(L$condt))[2]], 
                                               1, function(x) mean(x == 0)))
      fraczero$fraczerodiff <- abs(fraczero$fraczero1 - fraczero$fraczero2)
      vartpm <- data.frame(vartpm = apply(L$tpm, 1, var))
      # df <- merge(merge(merge(merge(melt(as.matrix(pvals)), avetpm, by.x = "Var1", by.y = 0, all = TRUE),
      #                         avecount, by.x = "Var1", by.y = 0, all = TRUE),
      #                   fraczero, by.x = "Var1", by.y = 0, all = TRUE),
      #             vartpm, by.x = "Var1", by.y = 0, all = TRUE)
      # df <- subset(df, rowSums(is.na(df)) == 0)
      # for (x in c("avetpm", "avecount", "fraczero", "vartpm", "fraczerodiff")) {
      #   print(ggplot(df, aes_string(x = x, y = "value", group = "Var2", col = "Var2")) + 
      #           geom_smooth() + theme_bw())
      # }
      df2 <- merge(merge(merge(merge(melt(as.matrix(padjs)), avetpm, by.x = "Var1", by.y = 0, all = TRUE),
                              avecount, by.x = "Var1", by.y = 0, all = TRUE),
                        fraczero, by.x = "Var1", by.y = 0, all = TRUE),
                  vartpm, by.x = "Var1", by.y = 0, all = TRUE)
      df2 <- subset(df2, rowSums(is.na(df2)) == 0)
      df2$sign <- df2$value <= 0.05
      # for (x in c("avetpm", "avecount", "fraczero", "vartpm", "fraczerodiff")) {
      #   print(ggplot(df2, aes_string(x = x, y = "value", group = "Var2", col = "Var2")) + 
      #           geom_smooth() + theme_bw())
      # }
      for (y in c("avetpm", "avecount", "vartpm")) {
        print(ggplot(df2, aes_string(x = "Var2", y = paste0("log2(", y, ")"), 
                                     fill = "Var2", dodge = "sign", alpha = "sign")) + 
                geom_boxplot() + theme_bw() + scale_fill_manual(values = colvec) + 
                scale_alpha_manual(values = c(0.2, 0.8)) + 
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
                guides(alpha = guide_legend(override.aes = 
                                              list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                   colour = NA))))
        print(ggplot(df2, aes_string(x = paste0("log2(", y, ")"), fill = "Var2", alpha = "sign")) + 
                geom_density() + theme_bw() + scale_fill_manual(values = colvec) + 
                scale_alpha_manual(values = c(0.2, 0.8)) + facet_wrap(~Var2) + 
                guides(alpha = guide_legend(override.aes = 
                                              list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                   colour = NA))))
      }
      for (y in c("fraczero", "fraczerodiff")) {
        print(ggplot(df2, aes_string(x = "Var2", y = y, fill = "Var2", dodge = "sign", alpha = "sign")) + 
                geom_boxplot() + theme_bw() + scale_fill_manual(values = colvec) + 
                scale_alpha_manual(values = c(0.2, 0.8)) + 
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
                guides(alpha = guide_legend(override.aes = 
                                              list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                   colour = NA))))
        print(ggplot(df2, aes_string(x = y, fill = "Var2", alpha = "sign")) + 
                geom_density() + theme_bw() + scale_fill_manual(values = colvec) + 
                scale_alpha_manual(values = c(0.2, 0.8)) + facet_wrap(~Var2) + 
                guides(alpha = guide_legend(override.aes = 
                                              list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                   colour = NA))))
      }
    }
  }
  
}

plot_results <- function(cobra, colvec) {
  cobraperf <- calculate_performance(cobra, aspects = "overlap", 
                                     type_venn = "adjp", thr_venn = 0.05)
  overlap(cobraperf) <- overlap(cobraperf)[, order(colnames(overlap(cobraperf)))]
  ol <- as.matrix(overlap(cobraperf))
  ol[is.na(ol)] <- 0
  
  ## --------------------- Number of detections ----------------------------- ##
  print(reshape2::melt(ol, varnames = c("gene", "method")) %>% 
          group_by(method) %>% 
          dplyr::summarise(nsignif = sum(value)) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate")) %>%
          dplyr::mutate(nbr_samples = factor(nbr_samples, levels = as.character(unique(sort(as.numeric(as.character(nbr_samples))))))) %>%
          ggplot(aes(x = method, y = nsignif, color = method, shape = nbr_samples)) + 
          geom_point(size = 5) + theme_bw() + xlab("") + ylab("Number of detections") + 
          scale_color_manual(values = colvec) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  ## ----------------------- Jaccard distance ------------------------------- ##
  jacc <- (t(ol) %*% ol)/(nrow(ol) - t(1 - ol) %*% (1 - ol))
  ## Set NaNs (caused by no significant genes) to 0
  w <- which(colSums(ol) == 0)
  for (i in w) {
    for (j in w) {
      if (is.na(jacc[i, j])) jacc[i, j] <- 1
    }
  }
  
  mdsjacc <- cmdscale(1 - jacc, k = 2, eig = TRUE)
  mdsjaccx <- data.frame(mdsjacc$points)
  colnames(mdsjaccx) <- c("MDS_Jaccard_1", "MDS_Jaccard_2")
  mdsjaccx$mth <- rownames(mdsjaccx)
  
  
  ## ---------------------- Spearman correlations --------------------------- ##
  tmpmt <- pval(cobra)
  if (!(all(colnames(iCOBRA::score(cobra)) %in% colnames(tmpmt)))) {
    sdn <- setdiff(colnames(iCOBRA::score(cobra)), colnames(tmpmt))
    tmpmt <- merge(tmpmt, -iCOBRA::score(cobra)[, sdn, drop = FALSE], by = 0, all = TRUE)
    rownames(tmpmt) <- tmpmt$Row.names
    tmpmt$Row.names <- NULL
  }
  if (!(all(colnames(padj(cobra)) %in% colnames(tmpmt)))) {
    sdn <- setdiff(colnames(padj(cobra)), colnames(tmpmt))
    tmpmt <- merge(tmpmt, padj(cobra)[, sdn, drop = FALSE], by = 0, all = TRUE)
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
    separate(mth, into = c("method", "nbr_samples", "replicate"))
  
  all_methods <- unique(tmpinfo$method)
  all_sizes <- as.numeric(unique(tmpinfo$nbr_samples))
  all_replicates <- as.numeric(unique(tmpinfo$replicate))
  max_nreps <- max(as.numeric(tmpinfo$replicate))
  
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
  df$nsamples <- factor(df$nsamples, levels = as.character(unique(sort(as.numeric(as.character(df$nsamples))))))
  ## Actual number
  print(
    df %>% 
      ggplot(aes(x = nmethods, y = value, 
                 group = interaction(nsamples, replicate), color = nsamples)) + 
      geom_line(size = 2) + theme_bw() + 
      scale_color_manual(values = colvec[c(1,4,8,7,2,6,5,10,3,9,11,12)], 
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
      scale_color_manual(values = colvec[c(1,4,8,7,2,6,5,10,3,9,11,12)], 
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
          scale_color_manual(values = colvec, name = "Method") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))  
  
  for (sz in all_sizes) {
    print(jaccm %>%
            dplyr::filter(nbr_samples1 == sz) %>%
            dplyr::filter(nbr_samples2 == sz) %>%
            ggplot(aes(x = method1, y = value, color = method1)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + ylab("Jaccard index") + 
            scale_color_manual(values = colvec, name = "Method") + 
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
          scale_color_manual(values = colvec) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))  
  
  for (sz in all_sizes) {
    print(spdistm %>%
            dplyr::filter(nbr_samples1 == sz) %>%
            dplyr::filter(nbr_samples2 == sz) %>%
            ggplot(aes(x = method1, y = value, color = method1)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + ylab("Spearman correlation") + 
            scale_color_manual(values = colvec) + 
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
      for (j in 1:3) {
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
      keepmethods <- intersect(c(paste0(all_methods, ".", sz, ".1"),
                                 paste0(all_methods, ".", sz, ".2"),
                                 paste0(all_methods, ".", sz, ".3")),
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
      for (j in 1:3) {
        keepmethods <- intersect(paste0(all_methods, ".", sz, ".", j), 
                                 rownames(spdist))
        if (length(keepmethods) > 0) {
          hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
          plot(hcl)
          pheatmap(1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order], 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   scale = "none", main = "Spearman correlation", display_numbers = TRUE)
        }
      }
      ## All replicates
      keepmethods <- intersect(c(paste0(all_methods, ".", sz, ".1"),
                                 paste0(all_methods, ".", sz, ".2"),
                                 paste0(all_methods, ".", sz, ".3")),
                               rownames(spdist))
      if (length(keepmethods) > 0) {
        hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
        plot(hcl)
        pheatmap(1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order], 
                 cluster_rows = FALSE, cluster_cols = FALSE, 
                 scale = "none", main = "Spearman correlation", display_numbers = TRUE)
      }
    }
    
    ## UpSet plots
    for (sz in all_sizes) {
      for (j in 1:3) {
        cpl <- prepare_data_for_plot(cobraperf, keepmethods = paste0(all_methods, ".", sz, ".", j), 
                                     colorscheme = colvec, incloverall = FALSE)
        plot_upset_with_reordering(cpl, nintersects = 25)
      }
    }
  }
  if (length(all_sizes) > 1) {
    for (mth in all_methods) {
      cpl <- prepare_data_for_plot(cobraperf, keepmethods = c(paste0(mth, ".", all_sizes, ".1"),
                                                              paste0(mth, ".", all_sizes, ".2"),
                                                              paste0(mth, ".", all_sizes, ".3")), 
                                   colorscheme = colvec, incloverall = FALSE)
      plot_upset_with_reordering(cpl, nintersects = 25)
    }
  }
  if (max_nreps > 1) {
    for (mth in all_methods) {
      for (sz in all_sizes) {
        tmpmth <- intersect(paste0(mth, ".", sz, ".", 1:max_nreps), basemethods(cobraperf))
        if (length(tmpmth) > 1) {
          cpl <- prepare_data_for_plot(cobraperf, keepmethods = tmpmth, colorscheme = colvec, 
                                       incloverall = FALSE)
          plot_upset_with_reordering(cpl, nintersects = 25)
        }
      }
    }
  }
}


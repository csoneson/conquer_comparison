source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

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

plot_consistency <- function(cobra, colvec, summary_data = list()) {
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
    # for (sz in all_sizes) {
    #   for (j in all_replicates) {
    #     ## A specific dataset (number of samples/replicate combination)
    #     keepmethods <- intersect(paste0(all_methods, ".", sz, ".", j),
    #                              rownames(spdist))
    #     if (length(keepmethods) > 0) {
    #       mdssub <- cmdscale(spdist[keepmethods, keepmethods], k = 2, eig = TRUE)
    #       mdssubx <- data.frame(mdssub$points)
    #       colnames(mdssubx) <- c("MDS_Spearman_1", "MDS_Spearman_2")
    #       mdssubx$mth <- rownames(mdssubx)
    #       print(mdssubx %>% 
    #               separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
    #               ggplot(aes(x = MDS_Spearman_1, y = MDS_Spearman_2, 
    #                          col = method, pch = nbr_samples)) + 
    #               geom_point(size = 10) + theme_bw() + 
    #               scale_color_manual(values = colvec, guide = FALSE) + 
    #               geom_label_repel(aes(label = method), size = 7) + 
    #               theme(legend.position = "bottom"))
    #     }
    #   }
    #   
    #   ## All replicates, given number of samples
    #   keepmethods <- intersect(unlist(lapply(paste0(all_methods, ".", sz), 
    #                                          function(m) paste0(m, ".", all_replicates))), 
    #                            rownames(spdist))
    #   if (length(keepmethods) > 0) {
    #     mdssub <- cmdscale(spdist[keepmethods, keepmethods], k = 2, eig = TRUE)
    #     mdssubx <- data.frame(mdssub$points)
    #     colnames(mdssubx) <- c("MDS_Spearman_1", "MDS_Spearman_2")
    #     mdssubx$mth <- rownames(mdssubx)
    #     print(mdssubx %>% 
    #             separate(mth, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>% 
    #             ggplot(aes(x = MDS_Spearman_1, y = MDS_Spearman_2, 
    #                        col = method, pch = nbr_samples)) + 
    #             geom_point(size = 5) + theme_bw() + 
    #             scale_color_manual(values = colvec, guide = FALSE) + 
    #             geom_label_repel(aes(label = method)) + 
    #             theme(legend.position = "bottom"))
    #   }
    # }
    
    ## Hierarchical clustering and heatmap
    # for (sz in all_sizes) {
    #   for (j in all_replicates) {
    #     keepmethods <- intersect(paste0(all_methods, ".", sz, ".", j), 
    #                              rownames(spdist))
    #     if (length(keepmethods) > 0) {
    #       hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
    #       plot(hcl)
    #       hclspd <- 1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order]
    #       pheatmap(hclspd, 
    #                cluster_rows = FALSE, cluster_cols = FALSE, 
    #                scale = "none", main = "Spearman correlation", display_numbers = TRUE,
    #                annotation_col = data.frame(method = colnames(hclspd), 
    #                                            row.names = colnames(hclspd)), 
    #                annotation_colors = list(method = structure(colvec, names = paste0(names(colvec),
    #                                                                                   ".", sz, ".", j))),
    #                annotation_legend = FALSE, annotation_names_col = FALSE,
    #                annotation_row = data.frame(method = rownames(hclspd), 
    #                                            row.names = rownames(hclspd)),
    #                annotation_names_row = FALSE)
    #     }
    #   }
    #   ## All replicates
    #   keepmethods <- intersect(unlist(lapply(paste0(all_methods, ".", sz), 
    #                                          function(m) paste0(m, ".", all_replicates))), 
    #                            rownames(spdist))
    #   if (length(keepmethods) > 0) {
    #     hcl <- hclust(d = as.dist(spdist[keepmethods, keepmethods]), method = "complete")
    #     plot(hcl)
    #     hclspd <- 1 - spdist[keepmethods, keepmethods][hcl$order, hcl$order]
    #     pheatmap(hclspd, 
    #              cluster_rows = FALSE, cluster_cols = FALSE, 
    #              scale = "none", main = "Spearman correlation", display_numbers = TRUE,
    #              annotation_col = data.frame(method = colnames(hclspd), 
    #                                          row.names = colnames(hclspd)), 
    #              annotation_colors = 
    #                list(method = structure(rep(colvec, each = length(all_replicates)), 
    #                                        names = unlist(lapply(paste0(names(colvec), ".", sz), 
    #                                                              function(m) paste0(m, ".", 
    #                                                                                 all_replicates))))),
    #              annotation_legend = FALSE, annotation_names_col = FALSE,
    #              annotation_row = data.frame(method = rownames(hclspd), 
    #                                          row.names = rownames(hclspd)),
    #              annotation_names_row = FALSE)
    #   }
    # }
    
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
    # for (sz in all_sizes) {
    #   for (j in all_replicates) {
    #     c2 <- colvec
    #     names(c2) <- paste0(names(c2), ".", sz, ".", j)
    #     km <- paste0(all_methods, ".", sz, ".", j)
    #     cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
    #                                  colorscheme = c2[km], incloverall = FALSE)
    #     if (ncol(overlap(cpl)) > 0) {
    #       overlap_table <- overlap(cpl)
    #       m <- min(which(colSums(overlap_table) > 0))
    #       if (is.finite(m)) 
    #         overlap_table <- overlap_table[, c(m, setdiff(1:ncol(overlap_table), m))]
    #       m <- max(which(colSums(overlap_table) > 0))
    #       if (is.finite(m))
    #         overlap_table <- overlap_table[, c(setdiff(1:ncol(overlap_table), m), m)]
    #       plotorder <- colnames(overlap_table)[order(colSums(overlap_table), 
    #                                                  seq(1:ncol(overlap_table)),
    #                                                  decreasing = "true")]
    #       nsets <- ncol(overlap_table)
    #       nintersects <- 25
    #       sets.bar.color <- plotcolors(cpl)[plotorder]
    #       overlap_table$fraczero <- truth(cobratmp)[match(rownames(overlap_table), 
    #                                                       rownames(truth(cobratmp))), 
    #                                              paste0("fraczero.", sz, ".", j)]
    #       upset(overlap_table, nsets = nsets, nintersects = nintersects, 
    #             sets.bar.color = sets.bar.color, 
    #             order.by = "freq", decreasing = TRUE, boxplot.summary = "fraczero",
    #             set.metadata = list(data = data.frame(sets = km, 
    #                                                   mth = km,
    #                                                   row.names = km),
    #                                 plots = list(list(type = "matrix_rows", 
    #                                                   column = "mth", 
    #                                                   colors = c2[km], alpha = 0.25))))
    #     }
    #   }
    # }
    
    # for (sz in all_sizes) {
    #   for (j in all_replicates) {
    #     c2 <- colvec
    #     names(c2) <- paste0(names(c2), ".", sz, ".", j)
    #     km <- intersect(paste0(all_methods, ".", sz, ".", j), basemethods(cobraperf))
    #     cpl <- prepare_data_for_plot(cobraperf, keepmethods = km, 
    #                                  colorscheme = c2[km], incloverall = FALSE)
    #     if (ncol(overlap(cpl)) > 0) {
    #       for (k in all_methods) {
    #         tmp_cp <- as(cpl, "COBRAPerformance")
    #         tryCatch({
    #           overlap(tmp_cp) <- overlap(tmp_cp)[overlap(tmp_cp)[, paste0(k, ".", sz, ".", j)] == 0 & 
    #                                                rowSums(overlap(tmp_cp)) > 0, ]
    #           km <- intersect(paste0(setdiff(all_methods, k), ".", sz, ".", j), basemethods(cobraperf))
    #           cpl1 <- prepare_data_for_plot(tmp_cp, keepmethods = km, 
    #                                         colorscheme = c2[km], incloverall = FALSE)
    #           plot_upset_with_reordering(cpl1, nintersects = 25, 
    #                                      mainbar.y.label = paste0("Intersection Size (only genes ", 
    #                                                               "not called DE by ", k, ")"), 
    #                                      set.metadata = list(data = data.frame(sets = km, 
    #                                                                            mth = km,
    #                                                                            row.names = km),
    #                                                          plots = list(list(type = "matrix_rows", 
    #                                                                            column = "mth", 
    #                                                                            colors = c2[km], alpha = 0.25))))}, 
    #           error = function(e) message(paste0("No overlap data for ", k, ".", sz, ".", j)))
    #       }
    #     }
    #   }
    # }
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


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

plot_consistency <- function(cobra, colvec, exts, summary_data = list()) {
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  ## Concordance plots
  maxrank <- 1000
  
  ## All sample sizes, all replicates
  pconc <- pval(cobratmp)
  for (i in colnames(pconc)) {
    pconc[, i] <- order(pconc[, i])
  }
  allm <- unique(get_method(colnames(pconc)))
  concvals <- do.call(rbind, lapply(allm, function(mth) {
    tmp <- pconc[, which(get_method(colnames(pconc)) == mth)]
    concval <- data.frame(t(sapply(1:maxrank, function(i) {
      p1 <- sum(table(unlist(tmp[1:i, ])) == ncol(tmp))
      p0.5 <- sum(table(unlist(tmp[1:i, ])) >= 0.5*ncol(tmp))
      c(k = i, p1 = p1, p0.5 = p0.5)
    })), stringsAsFactors = FALSE)
    concval$method <- mth
    concval
  }))
  print(ggplot(concvals, aes(x = k, y = p1, group = method, color = method)) + 
          geom_line() + scale_color_manual(values = colvec) + 
          theme_bw())
  print(ggplot(concvals, aes(x = k, y = p0.5, group = method, color = method)) + 
          geom_line() + scale_color_manual(values = colvec) + 
          theme_bw())
  summary_data$concordance_fullds <- rbind(summary_data$concordance_fullds, concvals)
  
  conc_auc <- concvals %>% dplyr::group_by(method) %>% 
    dplyr::summarize(auc1 = caTools::trapz(c(k, k[length(k)]), c(p1, 0))/(maxrank^2/2), 
                     auc0.5 = caTools::trapz(c(k, k[length(k)]), c(p0.5, 0))/(maxrank^2/2))
  summary_data$concordance_fullds_auc <- rbind(summary_data$concordance_fullds_auc, conc_auc)
  
  ## Per sample size
  allss <- unique(get_nsamples(colnames(pconc)))
  concvals_ss <- do.call(rbind, lapply(allm, function(mth) {
    do.call(rbind, lapply(allss, function(i) {
      tmp <- pconc[, intersect(which(get_method(colnames(pconc)) == mth),
                               which(get_nsamples(colnames(pconc)) == i)), drop = FALSE]
      if (ncol(tmp) > 1) {
        concval <- data.frame(t(sapply(1:maxrank, function(i) {
          p1 <- sum(table(unlist(tmp[1:i, ])) == ncol(tmp))
          p0.5 <- sum(table(unlist(tmp[1:i, ])) >= 0.5*ncol(tmp))
          c(k = i, p1 = p1, p0.5 = p0.5)
        })), stringsAsFactors = FALSE)
        concval$method <- mth
        concval$ncells <- i
        concval
      } else {
        NULL
      }
    }))
  }))
  concvals_ss$ncells <- factor(concvals_ss$ncells,
                               levels = sort(unique(as.numeric(as.character(concvals_ss$ncells)))))
  print(ggplot(concvals_ss, aes(x = k, y = p1, group = method, color = method)) + 
          geom_line() + scale_color_manual(values = colvec, name = "") + 
          theme_bw() + facet_wrap(~ncells) + theme(legend.position = "bottom"))
  print(ggplot(concvals_ss, aes(x = k, y = p0.5, group = method, color = method)) + 
          geom_line() + scale_color_manual(values = colvec, name = "") + 
          theme_bw() + facet_wrap(~ncells) + theme(legend.position = "bottom"))
  summary_data$concordance_byncells <- rbind(summary_data$concordance_byncells, concvals_ss)
  
  conc_auc_ss <- concvals_ss %>% dplyr::group_by(method, ncells) %>% 
    dplyr::summarize(auc1 = caTools::trapz(c(k, k[length(k)]), c(p1, 0))/(maxrank^2/2), 
                     auc0.5 = caTools::trapz(c(k, k[length(k)]), c(p0.5, 0))/(maxrank^2/2))
  summary_data$concordance_byncells_auc <- 
    rbind(summary_data$concordance_byncells_auc, conc_auc_ss)
  
  ## Pairwise, within sample size
  concvals_pairwise <- do.call(rbind, lapply(allm, function(mth) {
    do.call(rbind, lapply(allss, function(i) {
      tmp <- pconc[, intersect(which(get_method(colnames(pconc)) == mth),
                               which(get_nsamples(colnames(pconc)) == i)), drop = FALSE]
      if (ncol(tmp) > 1) {
        concval <- NULL
        for (j1 in 1:(ncol(tmp) - 1)) {
          for (j2 in (j1 + 1):(ncol(tmp))) {
            cv <- data.frame(t(sapply(1:maxrank, function(i) {
              p1 <- sum(table(unlist(tmp[1:i, c(j1, j2)])) == 2)
              c(k = i, p1 = p1)
            })), stringsAsFactors = FALSE)
            cv$ncells1 <- i
            cv$ncells2 <- i
            cv$replicate1 <- get_repl(colnames(tmp))[j1]
            cv$replicate2 <- get_repl(colnames(tmp))[j2]
            concval <- rbind(concval, cv)
          }
        }
        concval$method <- mth
        concval
      } else {
        NULL
      }
    }))
  }))
  summary_data$concordance_pairwise <- rbind(summary_data$concordance_pairwise, concvals_pairwise)
  
  conc_auc_pw <- concvals_pairwise %>% 
    dplyr::group_by(method, ncells1, ncells2, replicate1, replicate2) %>% 
    dplyr::summarize(auc1 = caTools::trapz(c(k, k[length(k)]), c(p1, 0))/(maxrank^2/2))
  summary_data$concordance_pairwise_auc <- 
    rbind(summary_data$concordance_pairwise_auc, conc_auc_pw)
  
  ## Between pairs of methods, within sample size/replicate
  allrepl <- unique(get_repl(colnames(pconc)))
  concvals_btwmth <- do.call(rbind, lapply(allss, function(ss) {
    do.call(rbind, lapply(allrepl, function(i) {
      tmp <- pconc[, intersect(which(get_repl(colnames(pconc)) == i),
                               which(get_nsamples(colnames(pconc)) == ss)), drop = FALSE]
      if (ncol(tmp) > 1) {
        concval <- NULL
        for (j1 in 1:(ncol(tmp) - 1)) {
          for (j2 in (j1 + 1):(ncol(tmp))) {
            cv <- data.frame(t(sapply(1:maxrank, function(i) {
              p1 <- sum(table(unlist(tmp[1:i, c(j1, j2)])) == 2)
              c(k = i, p1 = p1)
            })), stringsAsFactors = FALSE)
            cv$method1 <- get_method(colnames(tmp))[j1]
            cv$method2 <- get_method(colnames(tmp))[j2]
            concval <- rbind(concval, cv)
          }
        }
        concval$ncells <- ss
        concval$repl <- i
        concval
      } else {
        NULL
      }
    }))
  }))
  summary_data$concordance_betweenmethods <- 
    rbind(summary_data$concordance_betweenmethods, concvals_btwmth)
  
  conc_auc_btwmth <- concvals_btwmth %>% 
    dplyr::group_by(method1, method2, ncells, repl) %>% 
    dplyr::summarize(auc1 = caTools::trapz(c(k, k[length(k)]), c(p1, 0))/(maxrank^2/2))
  summary_data$concordance_betweenmethods_auc <- 
    rbind(summary_data$concordance_betweenmethods_auc, conc_auc_btwmth)
  
  ## Visualize cross-method consistency (average pairwise AUC)
  cmcons <- conc_auc_btwmth %>% dplyr::group_by(method1, method2) %>%
    dplyr::summarize(auc1 = mean(auc1)) %>% as.data.frame()
  cmcons <- rbind(cmcons, data.frame(method1 = unique(c(cmcons$method1, cmcons$method2)),
                                     method2 = unique(c(cmcons$method1, cmcons$method2)),
                                     auc1 = 1, stringsAsFactors = FALSE))
  cmcons <- dcast(cmcons, method1 ~ method2, value.var = "auc1")
  rownames(cmcons) <- cmcons$method1
  cmcons$method1 <- NULL
  stopifnot(all(rownames(cmcons)==colnames(cmcons)))
  for (i in 1:nrow(cmcons)) {
    for (j in 1:ncol(cmcons)) {
      if (is.na(cmcons[i, j])) cmcons[i, j] <- cmcons[j, i]
    }
  }
  pheatmap(cmcons, cluster_rows = TRUE, cluster_cols = TRUE)
  
  ## Visualize distributions of cross-method consistency
  cmdist <- conc_auc_btwmth
  for (i in 1:nrow(cmdist)) {
    if (cmdist[i, "method1"] > cmdist[i, "method2"]) {
      cmdist[i, c("method1", "method2")] <- cmdist[i, c("method2", "method1")]
    }
  }
  ggplot(cmdist, aes(x = auc1)) + geom_line(stat = "density") + 
    facet_grid(method1 ~ method2, scales = "free") + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank()) + 
    xlim(0, 1) + xlab("AUC")
  ggplot(cmdist, aes(x = auc1)) + geom_histogram() + 
    facet_grid(method1 ~ method2, scales = "free") + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank()) + 
    xlim(0, 1) + xlab("AUC")
  
  ## Get order to visualize distributions in
  hcl <- hclust(as.dist(1 - cmcons))
  hclord <- hcl$labels[hcl$order]
  cmdist2 <- cmdist %>% dplyr::group_by(method1, method2) %>%
    dplyr::arrange(auc1) %>% dplyr::mutate(ordr = 1:length(auc1), ycoord = 1)
  for (i in 1:nrow(cmdist2)) {
    if (match(cmdist2[i, "method1"], hclord) > match(cmdist2[i, "method2"], hclord)) {
      cmdist2[i, c("method1", "method2")] <- cmdist2[i, c("method2", "method1")]
    }
  }
  cmdist2$method1 <- factor(cmdist2$method1, levels = hclord[hclord %in% cmdist2$method1])
  cmdist2$method2 <- factor(cmdist2$method2, levels = hclord[hclord %in% cmdist2$method2])
  ggplot(cmdist2, aes(x = ordr, y = ycoord)) + geom_raster(aes(fill = auc1)) + 
    facet_grid(method1 ~ method2) + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank(),
          panel.border = element_blank(), strip.background = element_rect(colour = "white")) + 
    xlab("") + ylab("") + 
    scale_fill_continuous(low = "black", high = "yellow", name = "AUC")
  
  cmdist2 <- dplyr::mutate(cmdist2, ncells = as.numeric(as.character(ncells)))
  ggplot(cmdist2, aes(x = ncells, y = auc1)) + geom_point(size = 0.5) + 
    geom_smooth() + 
    facet_grid(method1 ~ method2) + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank(),
          strip.background = element_rect(colour = "white")) + 
    xlab("Number of cells per group") + ylab("AUC")
  
  ## Overlaps
  cobraperf <- calculate_performance(cobratmp, aspects = "overlap", 
                                     type_venn = "adjp", thr_venn = 0.05)
  overlap(cobraperf) <- overlap(cobraperf)[, order(colnames(overlap(cobraperf)))]
  ol <- as.matrix(overlap(cobraperf))
  ol[is.na(ol)] <- 0
  
  all_methods <- unique(get_method(colnames(ol)))
  all_sizes <- unique(get_nsamples(colnames(ol)))
  all_replicates <- unique(get_repl(colnames(ol)))
  max_nreps <- max(as.numeric(all_replicates))
  if (length(all_methods) > 1) {
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


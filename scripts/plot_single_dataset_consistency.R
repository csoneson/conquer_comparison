source("scripts/help_function_crossmethod_concordance.R")

#' Make UpSet plot, making sure that the first and last columns have at least
#' one significant feature
#' 
plot_upset_with_reordering <- function(cobraplot, nintersects, ...) {
  ## Reorder so that the first and last columns have something significant, if
  ## possible
  m <- min(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m)) 
    overlap(cobraplot) <- overlap(cobraplot)[, c(m, setdiff(seq_len(ncol(overlap(cobraplot))), m)), drop = FALSE]
  m <- max(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m))
    overlap(cobraplot) <- overlap(cobraplot)[, c(setdiff(seq_len(ncol(overlap(cobraplot))), m), m), drop = FALSE]
  tryCatch(plot_upset(cobraplot, order.by = "freq", decreasing = TRUE, nintersects = nintersects, ...),
           error = function(e) NULL)
}

plot_consistency <- function(cobra, concordances, colvec, exts, summary_data = list()) {
  ## ------------------------- Concordance plots ---------------------------- ##
  print(ggplot(concordances$concordance_fullds_bymethod %>% 
                 dplyr::ungroup() %>%
                 dplyr::mutate(method = gsub(exts, "", method)), 
               aes(x = k, y = nbr_genes, group = method, color = method)) +
          geom_line() + 
          scale_color_manual(values = structure(
            colvec, names = gsub(exts, "", names(colvec))), name = "") +
          theme_bw() + xlab("Number of top-ranked genes") + 
          ylab("Number of genes shared between all instances") + 
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 13)))

  for (k0 in c(100, 1000)) {
    if (max(concordances$concordance_betweenmethods_pairwise$k) >= k0) {
      help_function_crossmethod_concordance(concordances$concordance_betweenmethods_pairwise %>%
                                              dplyr::ungroup() %>%
                                              dplyr::mutate(method1 = gsub(exts, "", method1),
                                                            method2 = gsub(exts, "", method2)),
                                            k0 = k0)
    }
  }
  
  ## ------------------------------ Overlaps -------------------------------- ##
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
  cobraperf <- calculate_performance(cobratmp, aspects = "overlap", 
                                     type_venn = "adjp", thr_venn = 0.05)
  overlap(cobraperf) <- overlap(cobraperf)[, order(colnames(overlap(cobraperf))), drop = FALSE]
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
        message(paste0(sz, ".", j))
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
      message(mth)
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
        message(paste0(mth, ".", sz))
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


source("scripts/help_function_crossmethod_concordance.R")

#' Make UpSet plot, making sure that the first and last columns have at least
#' one significant feature
#' 
plot_upset_with_reordering <- function(cobraplot, nintersects, ...) {
  ## Reorder so that the first and last columns have something significant, if
  ## possible
  m <- min(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m)) 
    overlap(cobraplot) <- overlap(cobraplot)[, c(m, setdiff(1:ncol(overlap(cobraplot)), m))]
  m <- max(which(colSums(overlap(cobraplot)) > 0))
  if (is.finite(m))
    overlap(cobraplot) <- overlap(cobraplot)[, c(setdiff(1:ncol(overlap(cobraplot)), m), m)]
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

  print(ggplot(concordances$concordance_byncells_bymethod %>%
                 dplyr::ungroup() %>% 
                 dplyr::mutate(method = gsub(exts, "", method)), 
               aes(x = k, y = nbr_genes, group = method, color = method)) +
          geom_line() + 
          scale_color_manual(values = structure(
            colvec, names = gsub(exts, "", names(colvec))), name = "") +
          theme_bw() + facet_wrap(~ncells) + theme(legend.position = "bottom") + 
          xlab("Number of top-ranked genes") + 
          ylab("Number of genes shared between all instances") + 
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 13)))

  for (k0 in c(100, 1000)) {
    if (max(concordances$concordance_betweenmethods_pairwise$k >= k0)) {
      help_function_crossmethod_concordance(concordances$concordance_betweenmethods_pairwise %>%
                                              dplyr::ungroup() %>%
                                              dplyr::mutate(method1 = gsub(exts, "", method1),
                                                            method2 = gsub(exts, "", method2)),
                                            k0 = k0)
    }
  }
  
  ## Plot fraction of top-k genes shared by each number of methods
  k0 <- 100
  print(concordances$nbrshared_betweenmethods_all %>% 
          dplyr::mutate(ncells = factor(ncells, levels = sort(unique(as.numeric(as.character(ncells)))))) %>%
          dplyr::filter(k == k0) %>% dplyr::group_by(ncells, repl) %>%
          dplyr::mutate(tot_nbr_unique_genes = sum(nbr_genes)) %>%
          dplyr::mutate(frac_genes = nbr_genes/tot_nbr_unique_genes) %>%
          dplyr::mutate(m = sapply(1:length(frac_genes), function(i) sum(frac_genes[i:length(frac_genes)]))) %>%
          dplyr::mutate(ncells_repl = interaction(ncells, repl)) %>%
          ggplot(aes(x = nbr_occ, y = m, group = ncells_repl, color = ncells)) + 
          geom_line() + theme_bw() + xlab("Number of DE methods (M)") + 
          ylab(paste0("Fraction of top ", k0, " genes shared by at least M DE methods")) + 
          scale_color_discrete(name = "Number of \ncells per group"))
  
  
  ## ------------------------------ Overlaps -------------------------------- ##
  cobratmp <- cobra
  pval(cobratmp)[is.na(pval(cobratmp))] <- 1
  padj(cobratmp)[is.na(padj(cobratmp))] <- 1
  
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


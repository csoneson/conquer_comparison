summarize_timing <- function(figdir, datasets, exts, dtpext, cols,
                             singledsfigdir, cobradir, concordancedir, 
                             dschardir, origvsmockdir, distrdir, plotmethods, 
                             dstypes, pch_ncells) {

  ## Define ggplot2 parameters to use for all plots
  gglayers <- list(
    geom_boxplot(outlier.size = -1),
    theme_bw(),
    xlab(""),
    scale_color_manual(values = cols),
    guides(color = FALSE),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  )
  gglayersp <- c(gglayers, 
                 list(geom_point(position = position_jitter(width = 0.2), size = 0.5),
                      ylab("Relative computational \ntime requirement")))
  
  ## Initiate list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_timing", dtpext, ".pdf"), width = 10, height = 7)
  
  ## ------------------------------------------------------------------------ ##
  ## Read all timing information and calculate relative timing within each data
  ## set instance
  timing <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/timing/", ds, e, 
                     "_timing_summary_data.rds"))$timing_full %>%
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods) %>%
        dplyr::select(-user, -self)
    }))
  })) %>% group_by(dataset, filt, ncells, repl) %>%
    dplyr::mutate(rel_timing = timing/max(timing)) %>% ungroup() %>%
    dplyr::mutate(ncells_fact = factor(ncells, levels = 
                                         unique(sort(as.numeric(as.character(ncells)))))) %>%
    dplyr::ungroup()
  
  ## Read number of genes for each data set
  ngenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, "_nbr_called.rds")) %>% 
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::rename(ngenes = nbr_tested) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods) %>%
        dplyr::select(-nbr_sign_adjp0.05, -nbr_called, -nbr_nonsign_adjp0.05, 
                      -nbr_NA)
    }))
  }))
  timing <- dplyr::full_join(timing, ngenes) %>%
    dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::mutate(ngenes_cat = Hmisc::cut2(ngenes, g = 10))

  ## Add colors to the data frame
  timing$plot_color <- cols[as.character(timing$method)]

  ## Calculate exponent for single variable (keeping the other fixed)
  calc_expn <- function(y, x) {
    if (length(x) > 3) lm(log(y) ~ log(x))$coef[2]
    else NA
  }
  calc_mult <- function(y, x) {
    if (length(x) > 3) exp(lm(log(y) ~ log(x))$coef[1])
    else NA
  }
  exponfun <- function(x, a, p) {
    a * (x^p)
  }
  
  ## Dependence on number of cells for each method
  timing_ncelldep <- dplyr::full_join(
    timing %>% dplyr::group_by(method, ngenes_cat) %>% 
      dplyr::summarise(expn = calc_expn(y = timing, x = ncells)) %>%
      dplyr::ungroup(), 
    timing %>% dplyr::group_by(method, ngenes_cat) %>% 
      dplyr::summarise(mult = calc_mult(y = timing, x = ncells)) %>%
      dplyr::ungroup()
  )
  timing_ncelldep_curve <- 
    do.call(rbind, lapply(seq_len(nrow(timing_ncelldep)), function(i) {
      data.frame(method = timing_ncelldep[i, "method"], 
                 ngenes_cat = timing_ncelldep[i, "ngenes_cat"], 
                 x = 6:400, 
                 stringsAsFactors = FALSE) %>%
        dplyr::mutate(y = as.numeric(timing_ncelldep[i, "mult"]) * 
                        (x ^ as.numeric(timing_ncelldep[i, "expn"])))
  }))
  
  ## Dependence on number of genes for each method
  timing_ngenedep <- timing %>% group_by(method, ncells) %>% 
    dplyr::summarise(expn = calc_expn(y = timing, x = ngenes)) %>% ungroup()
  
  ## ------------------------------------------------------------------------ ##
  ## Boxplots
  ## Relative time requirement, all methods, ordered by median relative time
  ## requirement
  plots[["rel_timing_boxplot_comb"]] <- 
    ggplot(timing %>% 
             dplyr::mutate(method = forcats::fct_reorder(method, rel_timing, 
                                                         fun = median, na.rm = TRUE,
                                                         .desc = TRUE)), 
           aes(x = method, y = rel_timing, color = method)) + gglayersp + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = log10(5),
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
      hjust = 1, angle = 90) + 
    geom_hline(yintercept = 1.3, linetype = "dashed")
  print(plots[["rel_timing_boxplot_comb"]])

  ## Relative time requirement, all methods, ordered by median relative time
  ## requirement, log scale
  print(plots[["rel_timing_boxplot_comb"]] + scale_y_log10())
  
  ## Relative time requirement, all methods, ordered by median relative time
  ## requirement, facetted by number of cells
  print(plots[["rel_timing_boxplot_comb"]] + facet_wrap(~ ncells_fact))

  ## Exponent for modeling timing by number of genes, all methods, ordered by
  ## median exponent
  plots[["timing_exponent_ngenes"]] <- 
    ggplot(timing_ngenedep %>%
             dplyr::mutate(method = forcats::fct_reorder(method, expn, 
                                                         fun = median, na.rm = TRUE, 
                                                         .desc = TRUE)), 
           aes(x = method, y = expn, color = method)) + 
    gglayers + geom_point(position = position_jitter(width = 0.2), size = 1) + 
    ylab("Exponent in power model of \ntiming vs number of genes")
  print(plots[["timing_exponent_ngenes"]])
  
  ## Exponent for modeling timing by number of genes, all methods, ordered by
  ## median exponent
  plots[["timing_exponent_ncells"]] <- 
    ggplot(timing_ncelldep %>%
             dplyr::mutate(method = forcats::fct_reorder(method, expn, 
                                                         fun = median, na.rm = TRUE,
                                                         .desc = TRUE)), 
           aes(x = method, y = expn, color = method)) + 
    gglayers + geom_point(position = position_jitter(width = 0.2), size = 1) + 
    ylab("Exponent in power model of \ntiming vs number of cells") + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 1.9,
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 3, vjust = 0.5,
      hjust = 1, angle = 90) + 
    geom_hline(yintercept = 1.7, linetype = "dashed")
  print(plots[["timing_exponent_ncells"]])
  
  ## Timing dependence on number of cells, per method
  plots[["timing_dependence_permethod"]] <- 
    ggplot(timing, aes(x = ncells, y = timing, group = ngenes_cat, col = ngenes_cat)) + 
    facet_wrap(~ method, scales = "free") + geom_point(size = 0.5, alpha = 0.8) + 
    geom_line(data = timing_ncelldep_curve, aes(x = x, y = y), size = 1.25, alpha = 0.5) + 
    theme(legend.position = "bottom",
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12)) + 
    scale_color_manual(values = c("#332288", "#88CCEE", "#44AA99", "#117733",
                                  "#999933", "#DDCC77", "#661100", "#CC6677",
                                  "#882255", "#AA4499")) + 
    guides(color = guide_legend(nrow = 2, title = "Number of genes")) + 
    xlab("Number of cells per group") + ylab("Computational time requirement")
  print(plots[["timing_dependence_permethod"]])
  
  dev.off()
  
  pdf(paste0(figdir, "/summary_timing", dtpext, "_3d.pdf"), width = 15, height = 15)
  
  ## 3d plots 
  rbPal <- colorRampPalette(c('red','blue'))
  par(mfrow = c(3, 3))
  for (m in unique(timing$method)) {
    pl3d <- scatterplot3d(subset(timing, method == m)$ncells, xlab = "Number of cells per group",  
                          subset(timing, method == m)$ngenes, ylab = "Number of genes", 
                          subset(timing, method == m)$timing, zlab = "Computational time requirement", 
                          type = "h", pch = 19, main = m,
                          color = rbPal(20)[as.numeric(cut(subset(timing, method == m)$timing, 
                                                           breaks = 20))])
    model  <- lm(subset(timing, method == m)$timing ~ 
                   subset(timing, method == m)$ncells + subset(timing, method == m)$ngenes)
    pl3d$plane3d(model)
  }

  dev.off()
  
  ## --------------------------- Final summary plot ------------------------- ##
  pdf(paste0(figdir, "/timing_final", dtpext, ".pdf"), width = 12, height = 6)
  p <- plot_grid(plots[["rel_timing_boxplot_comb"]] + scale_y_log10(), 
                 plots[["timing_exponent_ncells"]],
                 labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
  print(p)
  dev.off()
  
  pdf(paste0(figdir, "/timing_final_ncellsdep", dtpext, ".pdf"), width = 14, height = 9.5)
  print(plots[["timing_dependence_permethod"]])
  dev.off()
  
  list(timing = timing %>% dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, 
                                         ngenes, ngenes_cat, timing, rel_timing),
       timing_ncelldep = timing_ncelldep)
}
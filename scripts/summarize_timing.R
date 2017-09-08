summarize_timing <- function(figdir, datasets, exts, dtpext, cols,
                             singledsfigdir, cobradir, concordancedir, 
                             dschardir, origvsmockdir, plotmethods) {

  ## Generate list to hold all plots
  plots <- list()
  
  pdf(paste0(figdir, "/summary_timing", dtpext, ".pdf"), width = 10, height = 7)
  
  ## Read all timing information
  timing <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/timing/", ds, e, 
                     "_timing_summary_data.rds"))$timing_full %>%
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
    }))
  }))
  ## Calculate relative timing within each data set instance
  timing <- timing %>% group_by(dataset, filt, ncells, repl) %>%
    dplyr::mutate(rel_timing = timing/max(timing)) %>% ungroup()
  
  ## Read number of genes for each data set
  ## 2d heatmaps
  ngenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, 
                     "_nbr_called.rds")) %>% 
        dplyr::mutate(repl = as.numeric(as.character(repl)),
                      ncells = as.numeric(as.character(ncells))) %>%
        dplyr::rename(ngenes = nbr_tested) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
    }))
  }))
  timing <- dplyr::full_join(timing, ngenes)
  
  timing$ncells_fact <- factor(timing$ncells, levels = 
                                 unique(sort(as.numeric(as.character(timing$ncells)))))
  ## Set plot symbols for number of cells per group
  ncells <- levels(timing$ncells_fact)
  pch <- c(16, 17, 15, 3, 7, 8, 4, 6, 9, 10, 11, 12, 13, 14)[1:length(ncells)]
  names(pch) <- as.character(ncells)
  
  ## Define colors for plotting
  cols <- structure(cols, names = gsub(paste(exts, collapse = "|"), "", names(cols)))
  
  ## Add colors and plot characters to the data frame
  timing$plot_color <- cols[as.character(timing$method)]
  timing$plot_char <- pch[as.character(timing$ncells_fact)]
  
  ## Boxplots
  plots[["rel_timing_boxplot_comb"]] <- 
    ggplot(timing, aes(x = method, y = rel_timing, color = method)) + 
    geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
    theme_bw() + xlab("") + ylab("Relative computational \ntime requirement") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = pch) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 3, title = "Number of \ncells per group", 
                                override.aes = list(size = 1.5)))
  print(plots[["rel_timing_boxplot_comb"]])

  plots[["rel_timing_boxplot_sep"]] <- 
    ggplot(timing, aes(x = method, y = rel_timing, color = method)) + 
    geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
    theme_bw() + xlab("") + ylab("Relative computational \ntime requirement") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = pch) + 
    facet_wrap(~ncells_fact) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 3, title = "Number of \ncells per group", 
                                override.aes = list(size = 1.5)))
  print(plots[["rel_timing_boxplot_sep"]])

  tmp <- timing %>% dplyr::mutate(method = as.character(method)) %>% 
    dplyr::group_by(method) %>% 
    dplyr::mutate(rel_timing_median = median(rel_timing, na.rm = TRUE)) %>%
    dplyr::ungroup()
  tmp$method <- factor(tmp$method, levels = unique(tmp$method[order(tmp$rel_timing_median, 
                                                                    decreasing = TRUE)]))
  plots[["rel_timing_boxplot_comb_log"]] <- 
    ggplot(tmp, aes(x = method, y = rel_timing, color = method)) + 
    geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
    theme_bw() + xlab("") + ylab("Relative computational \ntime requirement") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = pch) + 
    scale_y_log10() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 3, title = "Number of \ncells per group", 
                                override.aes = list(size = 1.5)))
  print(plots[["rel_timing_boxplot_comb_log"]])

  plots[["rel_timing_boxplot_sep_log"]] <- 
    ggplot(timing, aes(x = method, y = rel_timing, color = method)) + 
    geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = ncells_fact)) + 
    theme_bw() + xlab("") + ylab("Relative computational \ntime requirement") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = pch) + 
    scale_y_log10() + facet_wrap(~ncells_fact) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) + 
    guides(color = guide_legend(ncol = 2, title = ""),
           shape = guide_legend(ncol = 3, title = "Number of \ncells per group", 
                                override.aes = list(size = 1.5)))
  print(plots[["rel_timing_boxplot_sep_log"]])

  ## Barplots
  plots[["rel_timing_barplot_comb"]] <- 
    timing %>% group_by(method) %>% 
    dplyr::summarize(mean = mean(rel_timing), sd = sd(rel_timing)) %>%
    ggplot(aes(x = method, y = mean, fill = method)) + 
    geom_errorbar(aes(ymin = min(mean)/2, ymax = mean + sd), width = 0.2) + 
    geom_bar(stat = "identity") + 
    theme_bw() + xlab("") + ylab("Relative computational \ntime requirement") + 
    scale_fill_manual(values = cols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(plots[["rel_timing_barplot_comb"]])

  ## Help function to fit power model
  timing$ngenes_cat <- Hmisc::cut2(timing$ngenes, g = 10)
  ## Calculate exponent for single variable (keeping the other fixed)
  calc_expn <- function(y, x) {
    if (length(x) > 3) {
      lm(log(y) ~ log(x))$coef[2]
    } else {
      NA
    }
  }
  calc_mult <- function(y, x) {
    if (length(x) > 3) {
      exp(lm(log(y) ~ log(x))$coef[1])
    } else {
      NA
    }
  }
  exponfun <- function(x, a, p) {
    a * (x^p)
  }
  
  plots[["timing_exponent_ngenes"]] <- 
    timing %>% group_by(method, ncells) %>% 
    dplyr::summarise(expn = calc_expn(timing, ngenes)) %>%
    ggplot(aes(x = method, y = expn, color = method)) + 
    geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 1) + 
    theme_bw() + xlab("") + 
    ylab("Exponent in power model of \ntiming vs number of genes") + 
    scale_color_manual(values = cols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  print(plots[["timing_exponent_ngenes"]])
  
  tmp <- timing %>% dplyr::group_by(method, ngenes_cat) %>% 
    dplyr::summarise(expn = calc_expn(timing, ncells)) %>%
    dplyr::ungroup() %>% dplyr::mutate(method = as.character(method)) %>% 
    dplyr::group_by(method) %>%
    dplyr::mutate(expn_median = median(expn, na.rm = TRUE))
  tmp$method = factor(tmp$method, levels = unique(tmp$method[order(tmp$expn_median, 
                                                                   decreasing = TRUE)]))
  plots[["timing_exponent_ncells"]] <- 
    ggplot(tmp, aes(x = method, y = expn, color = method)) + geom_boxplot(outlier.size = -1) + 
    geom_point(position = position_jitter(width = 0.2), size = 1) + 
    theme_bw() + xlab("") + 
    ylab("Exponent in power model of \ntiming vs number of cells") + 
    scale_color_manual(values = structure(cols, names = gsub(paste(exts, collapse = "|"), 
                                                             "", names(cols))), name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  print(plots[["timing_exponent_ncells"]])
  
  ## Plot dependence on number of cells for each method
  tmp2 <- timing %>% dplyr::group_by(method, ngenes_cat) %>% 
    dplyr::summarise(mult = calc_mult(timing, ncells)) %>%
    dplyr::ungroup() %>% dplyr::mutate(method = as.character(method)) %>% 
    dplyr::group_by(method)
  tmp2 <- dplyr::full_join(tmp, tmp2)
  tmp3 <- do.call(rbind, lapply(1:nrow(tmp2), function(i) {
    data.frame(method = tmp2[i, "method"], 
               ngenes_cat = tmp2[i, "ngenes_cat"], 
               x = 6:100, 
               stringsAsFactors = FALSE) %>%
      dplyr::mutate(y = as.numeric(tmp2[i, "mult"]) * (x ^ as.numeric(tmp2[i, "expn"])))
  }))
  plots[["timing_dependence_permethod"]] <- 
    ggplot(timing, aes(x = ncells, y = timing, group = ngenes_cat, col = ngenes_cat)) + 
    facet_wrap(~ method, scales = "free") + geom_point(size = 0.5, alpha = 0.5) + 
    geom_line(data = tmp3, aes(x = x, y = y), size = 1.25) + 
    theme(legend.position = "bottom", #c(0.8, 0.075),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12)) + 
    scale_color_manual(values = c("#332288", "#88CCEE", "#44AA99", "#117733",
                                  "#999933", "#DDCC77", "#661100", "#CC6677",
                                  "#882255", "#AA4499")) + 
    guides(color = guide_legend(nrow = 2, title = "Number of genes")) + 
    xlab("Number of cells per group") + 
    ylab("Computational time requirement")
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
  p <- plot_grid(plot_grid(plots$rel_timing_boxplot_comb_log + theme(legend.position = "none"), 
                           plots$timing_exponent_ncells + theme(legend.position = "none"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 plot_grid(get_legend(plots$rel_timing_boxplot_comb_log + 
                                        theme(legend.position = "bottom") + 
                                        guides(colour = FALSE,
                                               shape = 
                                                 guide_legend(nrow = 3,
                                                              title = "Number of cells per group",
                                                              override.aes = list(size = 1.5),
                                                              title.theme = element_text(size = 12,
                                                                                         angle = 0),
                                                              label.theme = element_text(size = 10,
                                                                                         angle = 0),
                                                              keywidth = 1, default.unit = "cm"))),
                           NULL,
                           rel_widths = c(1, 1), nrow = 1),
                 rel_heights = c(1.7, 0.2), ncol = 1)
  print(p)
  dev.off()
  
  pdf(paste0(figdir, "/timing_final_ncellsdep", dtpext, ".pdf"), width = 14, height = 9.5)
  print(plots[["timing_dependence_permethod"]])
  dev.off()
  
  plots[c("rel_timing_boxplot_comb_log", "timing_exponent_ncells")]
}
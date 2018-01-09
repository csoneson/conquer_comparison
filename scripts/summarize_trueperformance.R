aspmod <- function(x) {
  if (x == "FDR") "FDP at\nadj.p = 0.05 cutoff"
  else if (x == "TPR") "TPR at\nadj.p = 0.05 cutoff"
  else if (x == "AUROC") "Area under ROC curve"
}

summarize_trueperformance <- function(figdir, datasets, exts, dtpext, cols,
                                      singledsfigdir, cobradir, concordancedir, 
                                      dschardir, origvsmockdir, distrdir, 
                                      plotmethods, dstypes, pch_ncells) {
  
  gglayers0 <- list(
    theme_bw(),
    xlab(""),
    scale_color_manual(values = cols),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)))
  gglayers <- c(list(geom_boxplot(outlier.size = -1),
                     geom_point(position = position_jitter(width = 0.2), size = 0.5)),
                gglayers0)
  gglayersfdr <- c(list(geom_hline(yintercept = 0.05)), 
                   gglayers,
                   scale_y_sqrt(breaks = c(0.05, 0.5, 1), limits = c(-0.1, 1.3)))
  
  ## Initialize list to hold all plots
  plots <- list()
  
  ## Define categorization of methods into "low FDP", "in range", "high FDP"
  getcat <- function(FDP, thr) {
    fracover <- length(which(FDP > thr))/length(FDP)
    fracunder <- length(which(FDP < thr))/length(FDP)
    fracoverx3 <- length(which(FDP > 3*thr))/length(FDP)
    fracunderx3 <- length(which(FDP < thr/3))/length(FDP)
    if (fracover >= 0.75 | fracoverx3 >= 0.5) "high FDP"
    else if (fracunder >= 0.75 | fracunderx3 >= 0.5) "low FDP"
    else "in range"
  }
  
  ## ------------------------------------------------------------------------ ##
  ## Read all true FDR/TPR information
  fdrtpr <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/performance_realtruth/", ds, e, 
                     "_performance_realtruth_summary_data.rds"))$FDRTPR
    }))
  }))
  fdrtpr_ihw <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/performance_realtruth/", ds, e, 
                     "_performance_realtruth_summary_data.rds"))$FDRTPR_IHW
    }))
  }))
  ## Read all true AUROC information
  auroc <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/performance_realtruth/", ds, e, 
                     "_performance_realtruth_summary_data.rds"))$AUROC
    }))
  }))
  ## Read all information about number of detected genes, fracNA etc
  nbrgenes <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(cobradir, "/", ds, e, "_nbr_called.rds")) %>%
        dplyr::mutate(method = paste(method, ncells, repl, sep = ".")) %>%
        dplyr::mutate(fracNA = nbr_NA/nbr_tested) %>%
        dplyr::select(method, dataset, filt, fracNA)
    }))
  }))
  fdrtprauc <- fdrtpr %>% dplyr::filter(thr == "thr0.05") %>%
    dplyr::left_join(nbrgenes) %>% 
    dplyr::left_join(dstypes, by = "dataset") %>% 
    dplyr::full_join(auroc) %>%
    tidyr::separate(method, c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
    dplyr::filter(method %in% plotmethods) %>%
    dplyr::group_by(thr, method, filt) %>% 
    dplyr::mutate(fdrcontrol = getcat(FDR, as.numeric(gsub("^thr", "", thr[1]))))  %>%
    dplyr::mutate(fdrcontrol = factor(fdrcontrol, levels = c("high FDP", "in range", "low FDP"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plot_color = cols[as.character(method)]) %>%
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells)))) %>%
    dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, thr, TPR, FDR, AUROC, fracNA, fdrcontrol)
  fdrtpr_ihw <- fdrtpr_ihw %>% dplyr::filter(thr == "thr0.05") %>%
    dplyr::left_join(nbrgenes) %>% 
    dplyr::left_join(dstypes, by = "dataset") %>% 
    tidyr::separate(method, c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
    dplyr::filter(method %in% plotmethods) %>%
    dplyr::group_by(thr, method, filt) %>% 
    dplyr::mutate(fdrcontrol = getcat(FDR, as.numeric(gsub("^thr", "", thr[1]))))  %>%
    dplyr::mutate(fdrcontrol = factor(fdrcontrol, levels = c("high FDP", "in range", "low FDP"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plot_color = cols[as.character(method)]) %>%
    dplyr::mutate(ncells_fact = factor(ncells, levels = sort(unique(ncells)))) %>%
    dplyr::select(method, dataset, dtype, filt, ncells_fact, repl, thr, TPR, FDR, fracNA, fdrcontrol)

  ## ----------------------------- Heatmaps --------------------------------- ##
  pdf(paste0(figdir, "/summary_trueperformance", dtpext, "_1.pdf"),
      width = 15, height = 4 * length(datasets))
  
  for (f in unique(fdrtprauc$filt)) {
    for (asp in c("FDR", "TPR", "AUROC")) {
      y <- fdrtprauc %>% 
        dplyr::filter(filt == f) %>%
        dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", ncells_fact, ".", repl)) %>%
        dplyr::select_("method", "dataset", asp) %>%
        reshape2::dcast(dataset ~ method, value.var = asp) %>%
        tidyr::separate(dataset, c("ds", "filt", "ncells", "repl"), sep = "\\.", remove = FALSE) %>%
        dplyr::arrange(ds, as.numeric(as.character(ncells))) %>% 
        dplyr::select(-ds, -filt, -ncells, -repl)  %>% as.data.frame()
      rownames(y) <- y$dataset
      y$dataset <- NULL
      
      annotation_row = data.frame(id = rownames(y)) %>% 
        tidyr::separate(id, c("dataset", "filt", "ncells", "repl"), sep = "\\.", remove = FALSE) %>%
        dplyr::mutate(
          ncells = factor(ncells, 
                          levels = as.character(sort(unique(as.numeric(as.character(ncells)))))))
      rownames(annotation_row) <- annotation_row$id
      
      pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", 
               main = paste0(aspmod(asp), " ", f),
               display_numbers = TRUE, 
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               breaks = seq(0, 1, length.out = 101), 
               annotation_row = dplyr::select(annotation_row, ncells, dataset), 
               show_rownames = FALSE,
               annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
               annotation_colors = list(method = cols[colnames(y)]),
               annotation_names_col = FALSE)
    }
  }
  dev.off()
  
  ## ------------------------------- Performance ------------------------------ ##
  pdf(paste0(figdir, "/summary_trueperformance", dtpext, "_2.pdf"),
      width = 10, height = 7)
  lineval <- c(FDR = 1, TPR = 1.1, AUROC = 1.1)
  yticks <- list(FDR = c(0.05, 0.5, 1), TPR = c(0, 0.5, 1), AUROC = c(0, 0.5, 1))
  
  for (f in unique(fdrtprauc$filt)) {
    for (asp in c("FDR", "TPR", "AUROC")) {
      ## All methods, order by median FDP
      p <- ggplot(fdrtprauc %>% dplyr::filter(filt == f) %>%
                    dplyr::mutate(method = forcats::fct_reorder(
                      method, FDR, fun = median, na.rm = TRUE, .desc = TRUE)),
                  aes_string(x = "method", y = asp, color = "method")) + 
        ylab(paste0(aspmod(asp))) + ggtitle(f)
      if (asp == "FDR") p <- p + gglayersfdr
      else p <- p + gglayers
      plots[[paste0(asp, "_all_", f)]] <- p
      print(plots[[paste0(asp, "_all_", f)]] + guides(color = FALSE))
      
      ## With n indicated
      x0 <- fdrtprauc %>% dplyr::filter(filt == f) %>%
        dplyr::mutate(method = forcats::fct_reorder(
          method, FDR, fun = median, na.rm = TRUE, .desc = TRUE))
      p <- ggplot(x0, aes_string(x = "method", y = asp, color = "method")) + 
        ylab(paste0(aspmod(asp))) + ggtitle(f)
      if (asp == "FDR") p <- p + gglayersfdr + 
        stat_summary(fun.data = function(x) {
          return(data.frame(y = 1,
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90)
      else p <- p + gglayers + 
        stat_summary(fun.data = function(x) {
          return(data.frame(y = 1.1,
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90) 
      p <- p + geom_hline(yintercept = lineval[asp], linetype = "dashed")
      plots[[paste0(asp, "_all_", f, "_withN")]] <- p
      print(plots[[paste0(asp, "_all_", f, "_withN")]] + guides(color = FALSE))
      
      ## All methods, order by fraction NA
      x0 <- fdrtprauc %>% dplyr::filter(filt == f) %>%
        dplyr::mutate(method = forcats::fct_reorder(method, fracNA, 
                                                    fun = median, na.rm = TRUE,
                                                    .desc = TRUE))
      p <- ggplot(x0, aes_string(x = "method", y = asp, color = "method")) + 
        ylab(paste0(aspmod(asp))) + ggtitle(f)
      if (asp == "FDR") p <- p + gglayersfdr + 
        stat_summary(fun.data = function(x) {
          return(data.frame(y = 1,
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90)
      else p <- p + gglayers + 
        stat_summary(fun.data = function(x) {
          return(data.frame(y = 1.1,
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90) 
      p <- p + geom_hline(yintercept = lineval[asp], linetype = "dashed")
      plots[[paste0(asp, "_all_fracnaorder_", f)]] <- p
      print(plots[[paste0(asp, "_all_fracnaorder_", f)]] + guides(color = FALSE))
      
      ## All methods, order by median FDP, facet by FDR control
      print(plots[[paste0(asp, "_all_", f)]] + 
              facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
              guides(color = FALSE))
      
      ## All methods, order by median FDP, facet by FDR control and data set
      print(plots[[paste0(asp, "_all_", f)]] + 
              facet_grid(dataset ~ fdrcontrol, scales = "free_x", space = "free_x") + 
              guides(color = FALSE))
      
      ## With IHW for all methods, order by median FDP, facet by FDR control
      if (asp %in% c("FDR", "TPR")) {
        x0 <- fdrtpr_ihw %>% dplyr::filter(filt == f) %>%
          dplyr::mutate(method = forcats::fct_reorder(method, FDR, 
                                                      fun = median, na.rm = TRUE,
                                                      .desc = TRUE))
        p <- ggplot(x0, aes_string(x = "method", y = asp, color = "method")) + 
          ylab(paste0(aspmod(asp))) + ggtitle(paste0(f, " IHW"))
        if (asp == "FDR") p <- p + gglayersfdr + 
          stat_summary(fun.data = function(x) {
            return(data.frame(y = 1,
                              label = paste0("n=", sum(!is.na(x)))))}, 
            geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
            hjust = -0.2, angle = 90)
        else p <- p + gglayers + 
          stat_summary(fun.data = function(x) {
            return(data.frame(y = 1.1,
                              label = paste0("n=", sum(!is.na(x)))))}, 
            geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
            hjust = -0.2, angle = 90) 
        p <- p + geom_hline(yintercept = lineval[asp], linetype = "dashed")
        plots[[paste0(asp, "_all_", f, "_IHW")]] <- p
        print(plots[[paste0(asp, "_all_", f, "_IHW")]] + 
                facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + guides(color = FALSE))
      }      
      
      ## Split by number of cells
      plots[[paste0(asp, "_byncells_sep_", f)]] <- 
        ggplot(fdrtprauc %>% dplyr::filter(filt == f),
               aes_string(x = "ncells_fact", y = asp, color = "method", group = "method")) + 
        geom_point(alpha = 0.25) + geom_smooth(se = FALSE) + 
        facet_wrap(~ dataset, scales = "free_x") + ylab(paste0(aspmod(asp))) + 
        gglayers0 + guides(color = guide_legend(ncol = 2, title = "")) + 
        ggtitle(f) + xlab("Number of cells per group")
      print(plots[[paste0(asp, "_byncells_sep_", f)]])
    }
    
    ## FDR vs TPR
    p2 <- fdrtprauc %>% dplyr::filter(filt == f) %>% 
      dplyr::group_by(dataset, ncells_fact, method) %>% 
      dplyr::summarize(TPR = median(TPR), FDR = median(FDR)) %>% 
      dplyr::ungroup() %>%
      ggplot(aes(x = FDR, y = TPR, color = method, label = method)) + 
      geom_point(size = 2) + 
      ylab(paste0("TPR at adj.p = 0.05 cutoff")) +
      facet_wrap(~ dataset + ncells_fact) + gglayers0 + 
      xlab(paste0("FDP at adj.p = 0.05 cutoff")) + 
      guides(color = guide_legend(ncol = 2, title = "")) + 
      ggtitle(f)
    plots[[paste0("fdrtpr_bydsncells_sep_", f)]] <- p2
    print(plots[[paste0("fdrtpr_bydsncells_sep_", f)]])
  }
  
  dev.off()

  ## -------------------------- Final summary plots ------------------------- ##
  pdf(paste0(figdir, "/trueperformance_final_byfdrcontrol", dtpext, ".pdf"), width = 12, height = 12)
  p <- plot_grid(plot_grid(plots[[paste0("FDR_all__withN")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("Without filtering"), 
                           plots[[paste0("FDR_all_TPM_1_25p_withN")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 plot_grid(plots[[paste0("TPR_all__withN")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("Without filtering") + 
                             scale_y_continuous(breaks = yticks$TPR, 
                                                limits = c(-0.01, 1.3)), 
                           plots[[paste0("TPR_all_TPM_1_25p_withN")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering") + 
                             scale_y_continuous(breaks = yticks$TPR, 
                                                limits = c(-0.01, 1.3)),
                           labels = c("C", "D"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 plot_grid(plots[[paste0("AUROC_all__withN")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("Without filtering") + 
                             scale_y_continuous(breaks = yticks$AUROC, 
                                                limits = c(-0.01, 1.3)), 
                           plots[[paste0("AUROC_all_TPM_1_25p_withN")]] +
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") +
                             theme(legend.position = "none") + 
                             ggtitle("After filtering") + 
                             scale_y_continuous(breaks = yticks$AUROC, 
                                                limits = c(-0.01, 1.3)),
                           labels = c("E", "F"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 rel_heights = c(1.7, 1.7, 1.7), ncol = 1)
  print(p)
  dev.off()
  
  ## Ordered by fracNA
  pdf(paste0(figdir, "/trueperformance_final_byfdrcontrol_fracnaorder", dtpext, ".pdf"), width = 12, height = 8)
  p <- plot_grid(plot_grid(plots[[paste0("FDR_all_fracnaorder_")]] + 
                             theme(legend.position = "none") + 
                             ggtitle("Without filtering"), 
                           plots[[paste0("FDR_all_fracnaorder_TPM_1_25p")]] + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 plot_grid(plots[[paste0("TPR_all_fracnaorder_")]] + 
                             theme(legend.position = "none") + 
                             ggtitle("Without filtering") + ylim(-0.01, 1.3), 
                           plots[[paste0("TPR_all_fracnaorder_TPM_1_25p")]] + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering") + ylim(-0.01, 1.3),
                           labels = c("C", "D"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 rel_heights = c(1.7, 1.7), ncol = 1)
  print(p)
  dev.off()
  
  ## IHW
  pdf(paste0(figdir, "/trueperformance_final_byfdrcontrol_ihw", dtpext, ".pdf"), width = 12, height = 8)
  p <- plot_grid(plot_grid(plots[[paste0("FDR_all__IHW")]] + 
                             theme(legend.position = "none") + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             ggtitle("Without filtering"), 
                           plots[[paste0("FDR_all_TPM_1_25p_IHW")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering"),
                           labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 plot_grid(plots[[paste0("TPR_all__IHW")]] + 
                             theme(legend.position = "none") + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             ggtitle("Without filtering") + ylim(-0.01, 1.3), 
                           plots[[paste0("TPR_all_TPM_1_25p_IHW")]] + 
                             facet_grid(~ fdrcontrol, scales = "free_x", space = "free_x") + 
                             theme(legend.position = "none") + 
                             ggtitle("After filtering") + ylim(-0.01, 1.3),
                           labels = c("C", "D"), align = "h", rel_widths = c(1, 1), nrow = 1),
                 rel_heights = c(1.7, 1.7), ncol = 1)
  print(p)
  dev.off()
  
  for (asp in c("FDR", "TPR", "AUROC")) {
    pdf(paste0(figdir, "/true", asp, "_final_sepbyds", dtpext, ".pdf"), width = 14, height = 6)
    p <- plot_grid(plot_grid(plots[[paste0(asp, "_byncells_sep_")]] + theme(legend.position = "none") + 
                               ggtitle("Without filtering") + ylim(-0.01, 1), 
                             plots[[paste0(asp, "_byncells_sep_TPM_1_25p")]] + 
                               theme(legend.position = "none") + 
                               ggtitle("After filtering") + ylim(-0.01, 1),
                             labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
                   get_legend(plots[[paste0(asp, "_byncells_sep_")]] + 
                                theme(legend.position = "bottom") + 
                                guides(colour = 
                                         guide_legend(nrow = 4,
                                                      title = "",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 12,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 10,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   rel_heights = c(1.7, 0.3), ncol = 1)
    print(p)
    dev.off()
    
    pdf(paste0(figdir, "/true", asp, "_final_sepbydsbox", dtpext, ".pdf"), width = 12, height = 8)
    p <- plot_grid(plots[[paste0(asp, "_all_")]] + 
                     facet_grid(dataset ~ fdrcontrol, scales = "free_x", space = "free_x") + 
                     guides(color = FALSE) + 
                     theme(legend.position = "none") + 
                     ggtitle("Without filtering") + ylim(-0.01, 1), 
                   plots[[paste0(asp, "_all_TPM_1_25p")]] + 
                     facet_grid(dataset ~ fdrcontrol, scales = "free_x", space = "free_x") + 
                     guides(color = FALSE) + 
                     theme(legend.position = "none") + 
                     ggtitle("After filtering") + ylim(-0.01, 1),
                   labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1)
    print(p)
    dev.off()
  }
  
  if (dtpext == "_sim") {
    write.table(fdrtprauc, file = "export_results/Figure4.csv", 
                row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
    
    pdf(paste0(figdir, "/trueFDR_for_slides", dtpext, "_filt.pdf"), width = 8, height = 4.8)
    print(plots[["FDR_all_TPM_1_25p"]] + guides(color = FALSE) + 
            ylab("FDP at adj.p = 0.05 cutoff") + scale_y_sqrt() + 
            theme(legend.position = "none") + ggtitle(""))
    dev.off()
  }
  
  fdrtprauc
}

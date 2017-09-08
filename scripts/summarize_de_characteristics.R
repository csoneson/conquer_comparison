summarize_de_characteristics <- function(figdir, datasets, exts, dtpext, cols,
                                         singledsfigdir, cobradir, concordancedir, 
                                         dschardir, origvsmockdir, plotmethods) {
  plots <- list()
  
  charname <- c(cvtpm = "CV(TPM)", fraczero = "Fraction zeros", 
                log2_avetpm = "log2(average TPM)", log2_vartpm = "log2(variance(TPM))",
                cvcpm = "CV(CPM)", log2_avecpm = "log2(average CPM)",
                log2_varcpm = "log2(variance(CPM))")
  
  pdf(paste0(figdir, "/summary_de_characteristics", dtpext, ".pdf"), 
      width = 12, height = 7)
  
  charac <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/results_characterization/", ds, e, 
                     "_results_characterization_summary_data.rds"))$stats_charac
    }))
  }))
  
  ## Change "mock" to "null" in data set names
  charac$dataset <- gsub("mock", "null", charac$dataset)
  
  ## Define colors for plotting
  cols <- structure(cols, names = gsub(paste(exts, collapse = "|"), "", names(cols)))
  
  for (f in unique(charac$filt)) {
    for (stat in c("tstat", "snr")) {
      x <- charac %>% dplyr::filter(filt == f) %>%
        dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
        dplyr::filter_(paste0("is.finite(", stat, ")")) %>% 
        dplyr::filter(charac != "fraczerodiff") %>%
        dplyr::filter(charac != "fraczeroround") %>%
        dplyr::filter(charac != "log2_avecount") %>%
        dplyr::filter(charac != "log2_avetpm") %>%
        dplyr::filter(charac != "log2_vartpm") %>%
        dplyr::filter(charac != "cvtpm") %>%
        tidyr::separate(Var2, into = c("method", "ncells", "repl"), sep = "\\.") %>%
        dplyr::mutate(charac = charname[charac]) %>%
        dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
        dplyr::filter(method %in% plotmethods)
      
      ## Visualize summary statistics for each characteristic
      statname <- switch(stat,
                         tstat = "t-statistic comparing significant\nand non-significant genes",
                         snr = "signal-to-noise statistic comparing\nsignificant and non-significant genes",
                         mediandiff = "median difference between\nsignificant and non-significant genes")
      p <- x %>% 
        ggplot(aes_string(x = "method", y = stat, color = "method")) + 
        geom_hline(yintercept = 0) + 
        geom_boxplot(outlier.size = -1) + 
        geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = dataset)) + 
        theme_bw() + 
        facet_wrap(~charac, scales = "free_y") + xlab("") + ylab(statname) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        scale_color_manual(values = cols) + 
        ggtitle(f) +
        guides(color = guide_legend(ncol = 2, title = ""),
               shape = guide_legend(ncol = 2, title = ""))
      plots[[paste0(stat, "_bystat_", f)]] <- p
      print(p)
      
      p <- x %>% 
        ggplot(aes_string(x = "charac", y = stat, color = "method", shape = "dataset")) + 
        geom_hline(yintercept = 0) + geom_point() + theme_bw() + 
        facet_wrap(~method, scales = "fixed") + xlab("") + ylab(statname) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        scale_color_manual(values = cols) + 
        ggtitle(f) + 
        guides(color = guide_legend(ncol = 2, title = ""),
               shape = guide_legend(ncol = 2, title = ""))
      plots[[paste0(stat, "_bymethod_", f)]] <- p
      print(p)
    
      p <- x %>% dplyr::group_by(charac) %>% 
        dplyr::mutate_(stat = paste0("(", stat, "-mean(", stat, "))/sd(", stat, ")")) %>% 
        dplyr::ungroup() %>% as.data.frame() %>%
        ggplot(aes_string(x = "charac", y = "stat", color = "method", shape = "dataset")) + 
        geom_hline(yintercept = 0) + 
        geom_point() + theme_bw() + facet_wrap(~method, scales = "fixed") + 
        xlab("") + ylab(paste0(statname, ",\ncentered and scaled across all instances)")) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        scale_color_manual(values = structure(cols, names = gsub(paste(exts, collapse = "|"),
                                                                 "", names(cols)))) + 
        ggtitle(f) + 
        guides(color = guide_legend(ncol = 2, title = ""),
               shape = guide_legend(ncol = 2, title = ""))
      plots[[paste0(stat, "_bymethod_scaled_", f)]] <- p
      print(p)
    }
  }
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  for (f in unique(charac$filt)) {
    for (stat in c("tstat", "snr")) {
      pdf(paste0(figdir, "/de_characteristics_final", ifelse(f == "", f, paste0("_", f)),
                 dtpext, "_", stat, ".pdf"), width = 10, height = 8)
      p <- plots[[paste0(stat, "_bystat_", f)]] + 
        theme(legend.position = "bottom") + 
        guides(colour = FALSE,
               shape = guide_legend(nrow = 2,
                                    title = "",
                                    override.aes = list(size = 1.5),
                                    title.theme = element_text(size = 12, angle = 0),
                                    label.theme = element_text(size = 12, angle = 0),
                                    keywidth = 1, default.unit = "cm"))
      print(p)
      dev.off()
    }
  }
  
  plots[c(paste0("tstat_bystat_", unique(charac$filt)),
          paste0("snr_bystat_", unique(charac$filt)))]
}
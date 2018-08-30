summarize_de_characteristics <- function(figdir, datasets, exts, dtpext, cols,
                                         singledsfigdir, cobradir, concordancedir, 
                                         dschardir, origvsmockdir, distrdir, 
                                         plotmethods, dstypes, pch_ncells) {
  
  gglayers <- list(
    geom_hline(yintercept = 0), 
    theme_bw(),
    xlab(""),
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          legend.position = "none"),
    scale_color_manual(values = cols)
  )
  
  plots <- list()
  
  charname <- c(cvtpm = "CV(TPM)", fraczero = "Fraction zeros", 
                log2_avetpm = "log2(average TPM)", log2_vartpm = "log2(variance(TPM))",
                cvcpm = "CV(CPM)", log2_avecpm = "log2(average CPM)",
                log2_varcpm = "log2(variance(CPM))",
                asinh_nb_minus_zinb_aic = "arcsinh(AIC[NB]-AIC[ZINB])")
  
  pdf(paste0(figdir, "/summary_de_characteristics", dtpext, ".pdf"), 
      width = 12, height = 9)
  
  charac <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      readRDS(paste0(singledsfigdir, "/results_characterization/", ds, e, 
                     "_results_characterization_summary_data.rds"))$stats_charac
    }))
  })) %>% dplyr::left_join(dstypes, by = "dataset") %>%
    dplyr::mutate(dataset = gsub("mock", "null", dataset)) %>%
    dplyr::filter(charac %in% c("cvcpm", "fraczero", "log2_avecpm", 
                                "log2_varcpm", "asinh_nb_minus_zinb_aic")) %>%
    tidyr::separate(Var2, into = c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(charac = charname[charac]) %>%
    dplyr::mutate(method = gsub(paste(exts, collapse = "|"), "", method)) %>%
    dplyr::filter(method %in% plotmethods) %>%
    dplyr::mutate(plot_color = cols[as.character(method)]) %>%
    dplyr::select(method, dataset, dtype, filt, ncells, repl, charac, snr)
    
  for (f in unique(charac$filt)) {
    for (stat in c("snr")) {
      x <- charac %>% dplyr::filter(filt == f) %>%
        dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
        dplyr::filter_(paste0("is.finite(", stat, ")"))

      ## Visualize summary statistics for each characteristic
      statname <- switch(stat,
                         tstat = "t-statistic comparing significant\nand non-significant genes",
                         snr = "signal-to-noise statistic comparing\nsignificant and non-significant genes",
                         mediandiff = "median difference between\nsignificant and non-significant genes")
      p <- x %>% dplyr::filter(!(charac == "arcsinh(AIC[NB]-AIC[ZINB])")) %>% 
        ggplot(aes_string(x = "method", y = stat, color = "method")) + 
        gglayers + geom_boxplot(outlier.size = -1) + 
        geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
        facet_wrap(~ charac, scales = "free_y") + xlab("") + ylab(statname) + 
        ggtitle(f)
      plots[[paste0(stat, "_bystat_", f)]] <- p
      print(p)
      
      ## With n indicated
      x0 <- x %>% dplyr::filter(!(charac == "arcsinh(AIC[NB]-AIC[ZINB])")) %>%
        dplyr::filter(!is.na(stat))
      p <- x0 %>% ggplot(aes_string(x = "method", y = stat, color = "method")) + 
        gglayers + geom_boxplot(outlier.size = -1) + 
        geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
        facet_wrap(~ charac, scales = "fixed") + xlab("") + ylab(statname) + 
        ggtitle(f) + stat_summary(fun.data = function(x) {
          return(data.frame(y = max(x0[["snr"]]),
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90) + 
        expand_limits(y = c(min(x0[[stat]]), 1.25*max(x0[[stat]]))) + 
        geom_hline(yintercept = max(x0[[stat]]), linetype = "dashed")
      plots[[paste0(stat, "_bystat_", f, "_withN")]] <- p
      print(p)
      
      ## Just a subset of the characteristics
      plots[[paste0(stat, "_bystat_", f, "_sub")]] <- 
        x %>% dplyr::filter(charac %in% c("log2(average CPM)", "Fraction zeros")) %>% 
        ggplot(aes_string(x = "method", y = stat, color = "method")) + 
        gglayers + geom_boxplot(outlier.size = -1) + 
        geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
        facet_wrap(~ charac, scales = "free_y") + xlab("") + ylab(statname) + 
        ggtitle(f)
      
      if ("arcsinh(AIC[NB]-AIC[ZINB])" %in% x$charac) {
        print(x %>% dplyr::filter(charac == "arcsinh(AIC[NB]-AIC[ZINB])") %>% 
                ggplot(aes_string(x = "method", y = stat, color = "method")) + 
                gglayers + geom_boxplot(outlier.size = -1) + 
                geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
                facet_wrap(~ charac, scales = "free_y") + xlab("") + ylab(statname) + 
                ggtitle(f))
      }
      
      p <- x %>% dplyr::filter(!(charac == "arcsinh(AIC[NB]-AIC[ZINB])")) %>%
        ggplot(aes_string(x = "charac", y = stat, color = "method")) + 
        gglayers + geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +  
        facet_wrap(~ method, scales = "fixed") + 
        ylab(statname) + ggtitle(f)
      plots[[paste0(stat, "_bymethod_", f)]] <- p
      print(p)
      
      ## With n indicated
      x0 <- x %>% dplyr::filter(!(charac == "arcsinh(AIC[NB]-AIC[ZINB])")) %>%
        dplyr::filter(!is.na(stat))
      p <- x0 %>% ggplot(aes_string(x = "charac", y = stat, color = "method")) + 
        gglayers + geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +  
        facet_wrap(~ method, scales = "fixed") + ylab(statname) + ggtitle(f) + 
        stat_summary(fun.data = function(x) {
          return(data.frame(y = max(x0[["snr"]]),
                            label = paste0("n=", sum(!is.na(x)))))}, 
          geom = "text", alpha = 1, color = "black", size = 2, vjust = 0.5, 
          hjust = -0.2, angle = 90) + 
        expand_limits(y = 1.75*max(x0[[stat]])) + 
        geom_hline(yintercept = max(x0[[stat]]), linetype = "dashed")
      plots[[paste0(stat, "_bymethod_", f, "_withN")]] <- p
      print(p)
    }
  }
  dev.off()
  
  ## -------------------------- Final summary plots ------------------------- ##
  for (f in unique(charac$filt)) {
    for (stat in c("snr")) {
      pdf(paste0(figdir, "/de_characteristics_final", ifelse(f == "", f, paste0("_", f)),
                 dtpext, "_", stat, ".pdf"), width = 10, height = 7)
      print(plots[[paste0(stat, "_bystat_", f, "_withN")]])
      dev.off()
      
      pdf(paste0(figdir, "/de_characteristics_final", ifelse(f == "", f, paste0("_", f)),
                 dtpext, "_", stat, "_bymethod.pdf"), width = 10, height = 10)
      print(plots[[paste0(stat, "_bymethod_", f, "_withN")]])
      dev.off()
    }
  }
  
  if (dtpext == "_real") {
    write.table(charac %>% dplyr::filter(filt == "") %>%
                  dplyr::filter(!is.na(snr)) %>%
                  dplyr::filter(is.finite(snr)) %>%
                  dplyr::filter(!(charac == "arcsinh(AIC[NB]-AIC[ZINB])")),
                file = "export_results/Figure2.csv", row.names = FALSE,
                col.names = TRUE, sep = ",", quote = FALSE)

    pdf(paste0(figdir, "/de_characteristics_for_slides",
               dtpext, "_snr.pdf"), width = 13, height = 6)
    print(plots[["snr_bystat__sub"]] + theme(strip.text = element_text(size = 15)))
    dev.off()
  }
  
  charac
}
summarize_filtering <- function(figdir, datasets, exts, dtpext, cols,
                                singledsfigdir, cobradir, concordancedir, 
                                dschardir, origvsmockdir, plotmethods) {
  exts <- setdiff(exts, "")
  
  for (e in exts) {
    pdf(paste0(figdir, "/summary_filtering", e, dtpext, ".pdf"), width = 12, height = 7)
    
    summary_data_list_orig <- lapply(datasets, function(ds) {
      readRDS(paste0(cobradir, "/", ds, "_nbr_called.rds"))
    })
    summary_data_list_filt <- lapply(datasets, function(ds) {
      readRDS(paste0(cobradir, "/", ds, e, "_nbr_called.rds"))
    })
    L <- rbind(do.call(rbind, summary_data_list_orig),
               do.call(rbind, summary_data_list_filt)) %>%
      dplyr::mutate(method = gsub(e, "", method)) %>%
      dplyr::filter(method == method[1]) %>%
      dplyr::group_by(dataset, ncells, repl) %>%
      dplyr::summarize(retain = nbr_tested[filt == gsub("^_", "", e)]/nbr_tested[filt == ""]) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ncells = factor(ncells, levels = sort(unique(as.numeric(as.character(ncells))))))
    
    print(ggplot(L, aes(x = dataset, y = retain)) + geom_boxplot(outlier.size = -1) + 
            geom_point(position = position_jitter(width = 0.2), aes(color = ncells)) + 
            theme_bw() + xlab("") + ylab("Retained fraction of original set of genes after filtering") + 
            guides(color = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_text(size = 13)))
    
    ## Plot library size vs fraction zeros across all data sets
    summary_data_list <- lapply(datasets, function(ds) {
      readRDS(paste0(dschardir, "/", ds, "_dataset_characteristics_summary_data.rds"))
    })
    L <- do.call(
      rbind, 
      lapply(summary_data_list, function(x) x$char_cells_m %>% 
               dplyr::filter(mtype %in% c("libsize", "fraczeroround")) %>%
               dplyr::mutate(cell = paste0(dataset, ".", cell, ".", ncells, ".", repl)) %>%
               dplyr::select(cell, mtype, value) %>%
               tidyr::spread(key = mtype, value = value) %>%
               tidyr::separate(cell, into = c("dataset", "cell", "ncells", "repl"), sep = "\\."))) %>%
      dplyr::mutate(ncells = factor(ncells, 
                                    levels = paste0(sort(unique(
                                      as.numeric(as.character(gsub(" cells per group", 
                                                                   "", ncells))))), 
                                      " cells per group")))
    print(ggplot(L, aes(x = libsize, y = fraczeroround, color = dataset, shape = ncells)) + 
            geom_point() + theme_bw() + xlab("Library size") + ylab("Fraction zeros after rounding") + 
            scale_shape_manual(values = 1:length(unique(L$ncells)), name = "") + 
            scale_color_discrete(name = "") + 
            theme(axis.text = element_text(size = 12), 
                  axis.title = element_text(size = 13)))
    
    print(ggplot(L, aes(x = libsize, y = fraczeroround, color = dataset, shape = ncells)) + 
            geom_point() + theme_bw() + xlab("Library size") + ylab("Fraction zeros after rounding") + 
            scale_shape_manual(values = 1:length(unique(L$ncells)), name = "") + 
            scale_color_discrete(name = "") + 
            scale_x_log10() + 
            theme(axis.text = element_text(size = 12), 
                  axis.title = element_text(size = 13)))
    
    dev.off()
  }
}
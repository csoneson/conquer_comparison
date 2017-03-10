summarize_filtering <- function(figdir, datasets, exts, dtpext, cols = cols) {
  pdf(paste0(figdir, "/summary_filtering", exts, dtpext, ".pdf"), width = 14, height = 7)
  
  summary_data_list_orig <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, "_nbr_called.rds"))
  })
  summary_data_list_filt <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, "_nbr_called.rds"))
  })
  L <- rbind(do.call(rbind, summary_data_list_orig),
             do.call(rbind, summary_data_list_filt)) %>%
    dplyr::mutate(method = gsub(exts, "", method)) %>%
    dplyr::filter(method == method[1]) %>%
    dplyr::group_by(dataset, ncells, repl) %>%
    dplyr::summarize(retain = nbr_tested[filt == gsub("^_", "", exts)]/nbr_tested[filt == ""]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ncells = factor(ncells, levels = sort(unique(as.numeric(as.character(ncells))))))
  
  print(ggplot(L, aes(x = dataset, y = retain)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = ncells)) + 
          theme_bw() + xlab("") + ylab("Retained fraction after filtering") + 
          scale_shape_discrete(name = "Number of cells") + 
          guides(shape = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
}
summarize_fracNA <- function(figdir, datasets, exts, dtpext, cols = cols) {
  ## -------------------------- Fraction NAs ---------------------------------- ##
  pdf(paste0(figdir, "/summary_fracNA", exts, dtpext, ".pdf"), width = 14, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, "_summary_data.rds"))
  })
  tmp_data <- lapply(summary_data_list, function(L) {
    L$all_data %>% filter(measurement == "fraczero") %>%
      dplyr::group_by(dataset, method, ncells, repl) %>%
      dplyr::summarize(fracNA = length(intersect(which(is.na(padj)), 
                                                 which(tested == TRUE)))/length(which(tested == TRUE)))
  })
  tmp_data <- do.call(rbind, tmp_data) %>% ungroup() %>% 
    dplyr::mutate(ncells = factor(ncells, levels = sort(unique(as.numeric(as.character(ncells))))))
  
  ## Remove extension from method name
  tmp_data$method <- gsub(exts, "", tmp_data$method)
  
  print(ggplot(tmp_data, aes(x = method, y = fracNA, color = method)) + 
          geom_boxplot(outlier.size = -1) +
          geom_point(position = position_jitter(width = 0.2), aes(shape = ncells)) + 
          theme_bw() + xlab("") + ylab("Fraction of NA adjusted p-values") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          scale_shape_discrete(name = "Number of cells") + 
          facet_wrap(~dataset) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  dev.off()
  
  pdf(paste0(figdir, "/summary_fracNA_comb", exts, dtpext, ".pdf"), width = 10, height = 7)
  print(ggplot(tmp_data, aes(x = method, y = fracNA, color = method)) + 
          geom_boxplot(outlier.size = -1) +
          geom_point(position = position_jitter(width = 0.2), aes(shape = ncells)) + 
          theme_bw() + xlab("") + ylab("Fraction of NA adjusted p-values") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          scale_shape_discrete(name = "Number of cells") + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  dev.off()
}
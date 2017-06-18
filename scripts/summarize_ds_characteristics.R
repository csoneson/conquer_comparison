summarize_ds_characteristics <- function(figdir, datasets, exts, dtpext, cols,
                                         singledsfigdir, cobradir, concordancedir, 
                                         dschardir, origvsmockdir, plotmethods) {
  
  X <- do.call(rbind, lapply(datasets, function(ds) {
    keepgroups <- fromJSON(file = paste0("config/", ds, ".json"))$keepgroups
    readRDS(paste0(dschardir, "/", ds, "_dataset_characteristics_summary_data.rds"))$char_cells_m %>%
      dplyr::filter(ncells == paste0(max(as.numeric(as.character(gsub(" cells per group", "", ncells)))), 
                                     " cells per group")) %>%
      dplyr::filter(mtype %in% c("fraczero", "libsize", "silhouette")) %>%
      dplyr::mutate(condition = condition == keepgroups[1])
  }))
  
  pdf(paste0(figdir, "/ds_characteristics_final", dtpext, ".pdf"), width = 18, height = 6)
  
  p1 <- ggplot(X %>% dplyr::filter(mtype == "fraczero"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + 
    theme_bw() + xlab("") + ylab("Fraction zeros per cell") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  p2 <- ggplot(X %>% dplyr::filter(mtype == "libsize"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + scale_y_log10() + 
    theme_bw() + xlab("") + ylab("Library size per cell") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  p3 <- ggplot(X %>% dplyr::filter(mtype == "silhouette"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + 
    theme_bw() + xlab("") + ylab("Silhouette width per cell") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))

  print(plot_grid(p1, p2, p3, labels = c("A", "B", "C"), align = "h", 
                  rel_widths = c(1, 1, 1), nrow = 1))
  
  dev.off()
}



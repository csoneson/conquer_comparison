summarize_pca <- function(figdir, datasets, exts) {
  ## ---------------------------------- PCA ----------------------------------- ##
  pdf(paste0(figdir, "/summary_pca", exts, ".pdf"), width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/results_characterization/", ds, exts, 
                   "_results_characterization_summary_data.rds"))
  })
  ## PCA of significant gene characteristics
  lapply(c(list(datasets), as.list(datasets)), function(ds) {
    for (stat in c("tstat", "mediandiff")) {
      x <- lapply(summary_data_list[ds], function(m) {
        m$stats_charac %>% dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
          dplyr::filter_(paste0("is.finite(", stat, ")")) %>% 
          dplyr::mutate(Var2 = paste0(Var2, ".", dataset, ".", filt)) %>%
          dplyr::select_("Var2", stat, "charac")  %>%
          dplyr::filter(charac != "fraczerodiff")
      })
      x <- do.call(rbind, x) %>% dcast(charac ~ Var2, value.var = stat)
      rownames(x) <- x$charac
      x$charac <- NULL
      x[is.na(x)] <- 0
      
      for (scl in c(TRUE, FALSE)) {
        pca <- prcomp(t(x), scale. = scl)
        annot <- data.frame(id = colnames(x), stringsAsFactors = FALSE) %>%
          tidyr::separate(id, into = c("method", "n_samples", "repl", "dataset", "filt"), 
                          sep = "\\.", remove = FALSE)
        print(ggplot(merge(annot, pca$x[, 1:2], by.x = "id", by.y = 0, all = TRUE), 
                     aes(x = PC1, y = PC2, color = method, shape = dataset)) +
                geom_point(size = 3) + theme_bw() + 
                scale_color_manual(values = cols) + 
                ggtitle(paste0(stat, ".scale=", scl)) + 
                guides(color = guide_legend(ncol = 2, title = ""),
                       shape = guide_legend(ncol = 2, title = "")))
        print(ggplot(merge(annot, pca$x[, 1:2], by.x = "id", by.y = 0, all = TRUE) %>%
                       dplyr::group_by(method, dataset) %>% 
                       dplyr::summarize(PC1 = mean(PC1), PC2 = mean(PC2)), 
                     aes(x = PC1, y = PC2, color = method, shape = dataset)) +
                geom_point(size = 4) + theme_bw() + 
                scale_color_manual(values = cols) + 
                ggtitle(paste0(stat, ".scale=", scl)) + 
                guides(color = guide_legend(ncol = 2, title = ""),
                       shape = guide_legend(ncol = 2, title = "")))
        print(ggplot(data.frame(id = rownames(pca$rotation), pca$rotation[, 1:2]), 
                     aes(x = PC1, y = PC2, label = id)) + geom_point() + geom_text_repel() +
                geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                             arrow = arrow(length = unit(0.03, "npc")), linetype = "dashed") + 
                theme_bw() + 
                ggtitle(paste0(stat, ".scale=", scl)))
        print(ggbiplot(pca, scale = 0, groups = annot$method, ellipse = TRUE,
                       ellipse.prob = 0.68, alpha = 0, var.axes = TRUE) + 
                theme_bw() + scale_color_manual(values = cols) + 
                ggtitle(paste0(stat, ".scale=", scl)) + 
                geom_point(data = merge(annot, pca$x[, 1:2], by.x = "id", by.y = 0, all = TRUE) %>%
                             dplyr::group_by(method) %>% 
                             dplyr::summarize(PC1 = mean(PC1), PC2 = mean(PC2)),
                           aes(x = PC1, y = PC2, col = method)) + 
                guides(color = guide_legend(ncol = 2, title = ""),
                       shape = guide_legend(ncol = 2, title = "")) + 
                coord_equal(ratio = diff(range(pca$x[, 1]))/diff(range(pca$x[, 2]))))
        plot(pca, main = "")
      }
    }
  })
  dev.off()
}
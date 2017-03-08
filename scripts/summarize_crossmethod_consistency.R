summarize_crossmethod_consistency <- function(figdir, datasets, exts, dtpext, cols = cols) {
  pdf(paste0(figdir, "/summary_crossmethod_consistency", exts, dtpext, ".pdf"), width = 14, height = 10)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/consistency/", ds, exts, 
                   "_consistency_summary_data.rds"))
  })
  concordances <- do.call(rbind, lapply(summary_data_list, 
                                        function(x) as.data.frame(x$concordance_betweenmethods_auc)))
  concordances$method1 <- gsub(exts, "", concordances$method1)
  concordances$method2 <- gsub(exts, "", concordances$method2)
  
  ## Visualize cross-method consistency (average pairwise AUC)
  cmcons <- concordances %>% dplyr::group_by(method1, method2) %>%
    dplyr::summarize(auc1 = mean(auc1)) %>% as.data.frame()
  cmcons <- rbind(cmcons, data.frame(method1 = unique(c(cmcons$method1, cmcons$method2)),
                                     method2 = unique(c(cmcons$method1, cmcons$method2)),
                                     auc1 = 1, stringsAsFactors = FALSE))
  cmcons <- dcast(cmcons, method1 ~ method2, value.var = "auc1")
  rownames(cmcons) <- cmcons$method1
  cmcons$method1 <- NULL
  stopifnot(all(rownames(cmcons)==colnames(cmcons)))
  for (i in 1:nrow(cmcons)) {
    for (j in 1:ncol(cmcons)) {
      if (is.na(cmcons[i, j])) cmcons[i, j] <- cmcons[j, i]
    }
  }
  pheatmap(cmcons, cluster_rows = TRUE, cluster_cols = TRUE)
  
  ## Visualize distributions of cross-method consistency
  cmdist <- concordances
  for (i in 1:nrow(cmdist)) {
    if (cmdist[i, "method1"] > cmdist[i, "method2"]) {
      cmdist[i, c("method1", "method2")] <- cmdist[i, c("method2", "method1")]
    }
  }
  print(ggplot(cmdist, aes(x = auc1)) + geom_line(stat = "density") + 
          facet_grid(method1 ~ method2, scales = "free") + theme_bw() + 
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
                strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank()) + 
          xlim(0, 1) + xlab("AUC"))
  print(ggplot(cmdist, aes(x = auc1)) + geom_histogram() + 
          facet_grid(method1 ~ method2, scales = "free") + theme_bw() + 
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
                strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank()) + 
          xlim(0, 1) + xlab("AUC"))
  
  ## Get order to visualize distributions in
  hcl <- hclust(as.dist(1 - cmcons))
  hclord <- hcl$labels[hcl$order]
  cmdist2 <- cmdist %>% dplyr::group_by(method1, method2) %>%
    dplyr::arrange(auc1) %>% dplyr::mutate(ordr = 1:length(auc1), ycoord = 1)
  for (i in 1:nrow(cmdist2)) {
    if (match(cmdist2[i, "method1"], hclord) > match(cmdist2[i, "method2"], hclord)) {
      cmdist2[i, c("method1", "method2")] <- cmdist2[i, c("method2", "method1")]
    }
  }
  cmdist2$method1 <- factor(cmdist2$method1, levels = hclord[hclord %in% cmdist2$method1])
  cmdist2$method2 <- factor(cmdist2$method2, levels = hclord[hclord %in% cmdist2$method2])
  print(ggplot(cmdist2, aes(x = ordr, y = ycoord)) + geom_raster(aes(fill = auc1)) + 
          facet_grid(method1 ~ method2) + theme_bw() + 
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
                strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank(),
                panel.border = element_blank(), strip.background = element_rect(colour = "white")) + 
          xlab("") + ylab("") + 
          scale_fill_continuous(low = "black", high = "yellow", name = "AUC"))
  
  cmdist2 <- dplyr::mutate(cmdist2, ncells = as.numeric(as.character(ncells)))
  print(ggplot(cmdist2, aes(x = ncells, y = auc1)) + geom_point(size = 0.5) + 
          geom_smooth() + 
          facet_grid(method1 ~ method2) + theme_bw() + 
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 6, angle = 90),
                strip.text.y = element_text(size = 6, angle = 0), panel.grid = element_blank(),
                strip.background = element_rect(colour = "white")) + 
          xlab("Number of cells per group") + ylab("AUC"))
  dev.off()
}
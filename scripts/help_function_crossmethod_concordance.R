suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(pheatmap))

help_function_crossmethod_concordance <- function(concordance_betweenmethods_auc) {
  ## Visualize cross-method consistency (average pairwise AUC)
  cmcons <- concordance_betweenmethods_auc %>% dplyr::group_by(method1, method2) %>%
    dplyr::summarize(auc = mean(auc)) %>% as.data.frame()
  cmcons <- rbind(cmcons, data.frame(method1 = unique(c(cmcons$method1, cmcons$method2)),
                                     method2 = unique(c(cmcons$method1, cmcons$method2)),
                                     auc = 1, stringsAsFactors = FALSE))
  cmcons <- dcast(cmcons, method1 ~ method2, value.var = "auc")
  rownames(cmcons) <- cmcons$method1
  cmcons$method1 <- NULL
  stopifnot(all(rownames(cmcons)==colnames(cmcons)))
  for (i in 1:nrow(cmcons)) {
    for (j in 1:ncol(cmcons)) {
      if (is.na(cmcons[i, j])) cmcons[i, j] <- cmcons[j, i]
    }
  }
  pheatmap(cmcons, cluster_rows = TRUE, cluster_cols = TRUE,
           main = "Average pairwise AUC across all data set instances", 
           color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "RdYlBu")))(100),
           breaks = seq(0, 1, length.out = 101))
  
  ## Visualize distributions of cross-method consistency
  cmdist <- concordance_betweenmethods_auc
  for (i in 1:nrow(cmdist)) {
    if (cmdist[i, "method1"] > cmdist[i, "method2"]) {
      cmdist[i, c("method1", "method2")] <- cmdist[i, c("method2", "method1")]
    }
  }
  print(ggplot(cmdist, aes(x = auc)) + geom_line(stat = "density") +
          facet_grid(method1 ~ method2, scales = "free") + theme_bw() +
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
                strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank()) +
          xlim(0, 1) + xlab("AUC"))
  print(ggplot(cmdist, aes(x = auc)) + geom_histogram() +
          facet_grid(method1 ~ method2, scales = "free") + theme_bw() +
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
                strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank()) +
          xlim(0, 1) + xlab("AUC"))
  
  ## Get order to visualize distributions in
  hcl <- hclust(as.dist(1 - cmcons))
  hclord <- hcl$labels[hcl$order]
  cmdist2 <- cmdist %>% dplyr::group_by(method1, method2) %>%
    dplyr::arrange(auc) %>% dplyr::mutate(ordr = 1:length(auc), ycoord = 1)
  for (i in 1:nrow(cmdist2)) {
    if (match(cmdist2[i, "method1"], hclord) > match(cmdist2[i, "method2"], hclord)) {
      cmdist2[i, c("method1", "method2")] <- cmdist2[i, c("method2", "method1")]
    }
  }
  cmdist2$method1 <- factor(cmdist2$method1, levels = hclord[hclord %in% cmdist2$method1])
  cmdist2$method2 <- factor(cmdist2$method2, levels = hclord[hclord %in% cmdist2$method2])
  print(ggplot(cmdist2, aes(x = ordr, y = ycoord)) + geom_raster(aes(fill = auc)) +
          facet_grid(method1 ~ method2) + theme_bw() +
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
                strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank(),
                panel.border = element_blank(), strip.background = element_rect(colour = "white")) +
          xlab("") + ylab("") +
          scale_fill_continuous(low = "black", high = "yellow", name = "AUC"))
  
  print(ggplot(cmdist2, aes(x = ordr, y = ycoord)) + geom_raster(aes(fill = auc)) +
          facet_grid(method1 ~ method2) + theme_bw() +
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
                strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank(),
                panel.border = element_blank(), strip.background = element_rect(colour = "white")) +
          xlab("") + ylab("") +
          scale_fill_continuous(low = "black", high = "yellow", name = "AUC", limits = c(0, 1)))
  
  cmdist2 <- dplyr::mutate(cmdist2, ncells = as.numeric(as.character(ncells)))
  print(ggplot(cmdist2, aes(x = ncells, y = auc)) + geom_point(size = 0.5) +
          geom_smooth() +
          facet_grid(method1 ~ method2) + theme_bw() +
          theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
                strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank(),
                strip.background = element_rect(colour = "white")) +
          xlab("Number of cells per group") + ylab("AUC"))
}
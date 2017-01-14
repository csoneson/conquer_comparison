args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggrepel))

datasets <- strsplit(datasets, ",")[[1]]

print(datasets)
print(filt)

cols <- c("#488d00", "#6400a6", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#a3c6ff", "#e6a900", "#a996ff", "#401600",
          "#ff6d9b", "#017671", "cyan", "red", "blue", "orange")
if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                        "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                        "voomlimma", "Wilcoxon", "BPSC", "MASTcounts", "MASTcountsDetRate", 
                        "MASTtpm"), exts)


summary_data_list <- lapply(datasets, function(ds) {
  readRDS(paste0("figures/summary_data/", ds, exts, "_summary_data.rds"))
})

pdf(paste0("figures/summary_crossds/summary_pca", exts, ".pdf"))

## PCA of significant gene characteristics
for (stat in c("tstat", "mediandiff")) {
  x <- lapply(summary_data_list, function(m) {
    m$stats_charac %>% dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
      dplyr::mutate(Var2 = paste0(Var2, ".", dataset, ".", filt)) %>%
      dplyr::select_("Var2", stat, "charac") 
  })
  x <- do.call(rbind, x) %>% dcast(charac ~ Var2, value.var = stat)
  rownames(x) <- x$charac
  x$charac <- NULL
  
  for (scl in c(TRUE, FALSE)) {
    pca <- prcomp(t(x), scale. = scl)
    annot <- data.frame(id = colnames(x), stringsAsFactors = FALSE) %>%
      tidyr::separate(id, into = c("method", "n_samples", "repl", "dataset", "filt"), 
                      sep = "\\.", remove = FALSE)
    print(ggplot(merge(annot, pca$x[, 1:2], by.x = "id", by.y = 0, all = TRUE), 
                 aes(x = PC1, y = PC2, color = method, shape = dataset)) +
            geom_point(size = 3) + theme_bw() + 
            scale_color_manual(values = cols) + 
            ggtitle(paste0(stat, ".scale=", scl)))
    print(ggplot(data.frame(id = rownames(pca$rotation), pca$rotation[, 1:2]), 
                 aes(x = PC1, y = PC2, label = id)) + geom_point() + geom_text_repel() +
            geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                         arrow = arrow(length = unit(0.03, "npc")), linetype = "dashed") + 
            theme_bw() + 
            ggtitle(paste0(stat, ".scale=", scl)))
    plot(pca, main = "")
  }
}
dev.off()

pdf(paste0("figures/summary_crossds/summary_heatmaps", exts, ".pdf"),
    width = 10, height = 3 * length(datasets))

## Heatmap of true FPRs (fraction of nominal p-values below 0.05)
y <- lapply(summary_data_list, function(m) {
  m$fracpbelow0.05 %>% 
    tidyr::separate(method, c("method", "n_samples", "repl"), sep = "\\.") %>%
    #dplyr::arrange(as.numeric(as.character(n_samples))) %>%
    dplyr::mutate(dataset = paste0(dataset, ".", filt, ".", n_samples, ".", repl)) %>%
    dplyr::select(method, dataset, FPR)
})
y <- do.call(rbind, y) %>% dcast(dataset ~ method, value.var = "FPR") %>%
  tidyr::separate(dataset, c("ds", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
  dplyr::arrange(ds, as.numeric(as.character(n_samples))) %>% 
  dplyr::select(-ds, -filt, -n_samples, -repl)  %>% as.data.frame()
rownames(y) <- y$dataset
y$dataset <- NULL

annotation_row = data.frame(id = rownames(y)) %>% 
  tidyr::separate(id, c("dataset", "filt", "n_samples", "repl"), sep = "\\.", remove = FALSE) %>%
  dplyr::mutate(n_samples = factor(n_samples, 
                                   levels = as.character(sort(unique(as.numeric(as.character(n_samples)))))))
rownames(annotation_row) <- annotation_row$id

pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", main = "True FPR",
         display_numbers = TRUE, 
         annotation_row = dplyr::select(annotation_row, n_samples, dataset), show_rownames = FALSE,
         annotation_col = data.frame(method = colnames(y), row.names = colnames(y)),
         annotation_colors = list(method = structure(cols, names = names(cols))[colnames(y)]),
         annotation_names_col = FALSE)
dev.off()

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
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(RColorBrewer))

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datasets)
print(filt)

cols <- c("#488d00", "#6400a6", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#a3c6ff", "#e6a900", "#a996ff", "#401600",
          "#ff6d9b", "#017671", "cyan", "red", "blue", "orange",
          "#777777", "#7BAFDE", "#F6C141", "#90C987", "#1965B0",
          "#882E72", "#F7EE55")
if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                        "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                        "voomlimma", "Wilcoxon", "BPSC", "MASTcounts", "MASTcountsDetRate", 
                        "MASTtpm", "zingeRauto", "Seurat", "DESeq2census", "edgeRLRTcensus",
                        "DESeq2nofilt", "Seuratnofilt", "NODESnofilt"), exts)


summary_data_list <- lapply(datasets, function(ds) {
  readRDS(paste0("figures/summary_data/", ds, exts, "_summary_data_orig_vs_mock.rds"))
})
jaccm0.05 <- do.call(rbind, lapply(summary_data_list, function(x) x$jaccard_adjp0.05))
jaccm0.1 <- do.call(rbind, lapply(summary_data_list, function(x) x$jaccard_adjp0.1))
jaccm0.2 <- do.call(rbind, lapply(summary_data_list, function(x) x$jaccard_adjp0.2))
spearm <- do.call(rbind, lapply(summary_data_list, function(x) x$spearman))

ttest <- function(x, y) {
  tryCatch({
    t.test(y[x == "original"], y[x == "mock"], var.equal = FALSE)$stat
  }, error = function(e) NA)
}

## T-statistic of Spearman correlations and Jaccard coefficients between original and mock
nbr_keep <- unique(intersect(subset(jaccm, tp == "original")$nbr_samples1,
                             subset(jaccm, tp == "mock")$nbr_samples1))
spearm %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
  dplyr::filter(nbr_samples1 == nbr_samples2) %>%
  dplyr::group_by(dataset, nbr_samples1, method1) %>% 
  dplyr::summarize(tstat = ttest(tp, value)) %>%
  ggplot(aes(x = method1, y = tstat, col = method1, shape = dataset)) + 
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_bw() + xlab("") + ylab("t-statistic, Spearman correlation (original - mock)") + 
  scale_color_manual(values = cols, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

jaccm0.05 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
  dplyr::filter(nbr_samples1 == nbr_samples2) %>%
  dplyr::group_by(dataset, nbr_samples1, method1) %>% 
  dplyr::summarize(tstat = ttest(tp, value)) %>%
  ggplot(aes(x = method1, y = tstat, col = method1, shape = dataset)) + 
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.05 (original - mock)") + 
  scale_color_manual(values = cols, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

jaccm0.1 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
  dplyr::filter(nbr_samples1 == nbr_samples2) %>%
  dplyr::group_by(dataset, nbr_samples1, method1) %>% 
  dplyr::summarize(tstat = ttest(tp, value)) %>%
  ggplot(aes(x = method1, y = tstat, col = method1, shape = dataset)) + 
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.1 (original - mock)") + 
  scale_color_manual(values = cols, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

jaccm0.2 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
  dplyr::filter(nbr_samples1 == nbr_samples2) %>%
  dplyr::group_by(dataset, nbr_samples1, method1) %>% 
  dplyr::summarize(tstat = ttest(tp, value)) %>%
  ggplot(aes(x = method1, y = tstat, col = method1, shape = dataset)) + 
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.2 (original - mock)") + 
  scale_color_manual(values = cols, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

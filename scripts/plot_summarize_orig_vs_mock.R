args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datasets)
print(filt)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)


summary_data_list <- lapply(datasets, function(ds) {
  readRDS(paste0("figures/orig_vs_mock/", ds, exts, "_orig_vs_mock_summary_data.rds"))
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

pdf(paste0("figures/summary_crossds/summary_orig_vs_mock", exts, ".pdf"), width = 10, height = 7)

## T-statistic of Spearman correlations and Jaccard coefficients between original and mock
nbr_keep <- unique(intersect(subset(jaccm0.05, tp == "original")$nbr_samples1,
                             subset(jaccm0.05, tp == "mock")$nbr_samples1))
print(spearm %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
        dplyr::filter(nbr_samples1 == nbr_samples2) %>%
        dplyr::group_by(dataset, nbr_samples1, method1) %>% 
        dplyr::summarize(tstat = ttest(tp, value)) %>% 
        ggplot(aes(x = method1, y = tstat, col = method1)) + 
        geom_boxplot(outlier.size = 0) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + ylab("t-statistic, Spearman correlation (original - mock)") + 
        scale_color_manual(values = cols, name = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

print(jaccm0.05 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
        dplyr::filter(nbr_samples1 == nbr_samples2) %>%
        dplyr::group_by(dataset, nbr_samples1, method1) %>% 
        dplyr::summarize(tstat = ttest(tp, value)) %>%
        ggplot(aes(x = method1, y = tstat, col = method1)) + 
        geom_boxplot(outlier.size = 0) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.05 (original - mock)") + 
        scale_color_manual(values = cols, name = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
      
print(jaccm0.1 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
        dplyr::filter(nbr_samples1 == nbr_samples2) %>%
        dplyr::group_by(dataset, nbr_samples1, method1) %>% 
        dplyr::summarize(tstat = ttest(tp, value)) %>%
        ggplot(aes(x = method1, y = tstat, col = method1)) + 
        geom_boxplot(outlier.size = 0) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.1 (original - mock)") + 
        scale_color_manual(values = cols, name = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

print(jaccm0.2 %>% dplyr::filter(nbr_samples1 %in% nbr_keep & nbr_samples2 %in% nbr_keep) %>%
        dplyr::filter(nbr_samples1 == nbr_samples2) %>%
        dplyr::group_by(dataset, nbr_samples1, method1) %>% 
        dplyr::summarize(tstat = ttest(tp, value)) %>%
        ggplot(aes(x = method1, y = tstat, col = method1)) + 
        geom_boxplot(outlier.size = 0) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + ylab("t-statistic, Jaccard index, adjp = 0.2 (original - mock)") + 
        scale_color_manual(values = cols, name = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

dev.off()

saveRDS(NULL, paste0("figures/summary_crossds/summary_orig_vs_mock", exts, ".rds"))

sessionInfo()


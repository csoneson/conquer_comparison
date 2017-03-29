args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datasets)
print(filt)
print(dtpext)
print(origvsmockdir)
print(figdir)

source("scripts/plot_setup.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

ttest <- function(x, y) {
  tryCatch({
    t.test(y[x == "signal"], y[x == "mock"], var.equal = FALSE)$stat
  }, error = function(e) NA)
}

plots <- list()

summary_data_list <- lapply(datasets, function(ds) {
  readRDS(paste0(origvsmockdir, "/", ds, exts, "_orig_vs_mock_summary_data.rds"))
})
concs <- do.call(rbind, lapply(summary_data_list, function(x) x$concordances))

pdf(paste0(figdir, "/summary_orig_vs_mock", exts, dtpext, ".pdf"), 
    width = 10, height = 7)

nbr_keep <- unique(intersect(subset(concs, tp == "signal")$ncells,
                             subset(concs, tp == "mock")$ncells))
## Remove extension from method name
concs$method <- gsub(exts, "", concs$method)

concsum <- concs %>% as.data.frame() %>%
  dplyr::filter(ncells %in% nbr_keep) %>%
  dplyr::group_by(dataset, ncells, method) %>% 
  dplyr::summarize(tstat = ttest(tp, AUCs),
                   mediandiff = median(AUCs[tp == "signal"]) - median(AUCs[tp == "mock"]))

p <- concs %>% dplyr::filter(tp == "signal") %>%
  dplyr::filter(ncells %in% nbr_keep) %>%
  ggplot(aes(x = method, y = AUCs, col = method)) + 
  geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
  geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) + 
  theme_bw() + xlab("") + ylab("Area under concordance curve, signal data set") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["auc_signal_comb"]] <- p

p <- concs %>% dplyr::filter(tp == "mock") %>%
  dplyr::filter(ncells %in% nbr_keep) %>%
  ggplot(aes(x = method, y = AUCs, col = method)) + 
  geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
  geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) + 
  theme_bw() + xlab("") + ylab("Area under concordance curve, mock data set") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["auc_mock_comb"]] <- p

p <- concs %>% dplyr::filter(tp == "signal") %>%
  dplyr::filter(ncells %in% nbr_keep) %>%
  ggplot(aes(x = method, y = AUCs, col = method)) + 
  geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  facet_wrap(~dataset) + 
  theme_bw() + xlab("") + ylab("Area under concordance curve, signal data set") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["auc_signal_sep"]] <- p

p <- concs %>% dplyr::filter(tp == "mock") %>%
  dplyr::filter(ncells %in% nbr_keep) %>%
  ggplot(aes(x = method, y = AUCs, col = method)) + 
  geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  facet_wrap(~dataset) + 
  theme_bw() + xlab("") + ylab("Area under concordance curve, signal data set") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["auc_mock_sep"]] <- p

p <- concsum %>% 
  ggplot(aes(x = method, y = sign(tstat) * sqrt(abs(tstat)), col = method)) + 
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
  theme_bw() + xlab("") + ylab("sqrt(t-statistic, area under concordance curve (signal - mock))") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["tstat_auc"]] <- p

p <- concsum %>% 
  ggplot(aes(x = method, y = mediandiff, col = method)) + 
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
  theme_bw() + xlab("") + 
  ylab("difference between median area under concordance curve (signal - mock)") + 
  scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
  scale_shape_discrete(name = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13))
print(p)
plots[["mediandiff_auc"]] <- p

dev.off()

saveRDS(plots, paste0(figdir, "/summary_orig_vs_mock", exts, dtpext, "_plots.rds"))

sessionInfo()


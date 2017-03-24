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
    t.test(y[x == "original"], y[x == "mock"], var.equal = FALSE)$stat
  }, error = function(e) NA)
}

summary_data_list <- lapply(datasets, function(ds) {
  readRDS(paste0(origvsmockdir, "/", ds, exts, "_orig_vs_mock_summary_data.rds"))
})
concs <- do.call(rbind, lapply(summary_data_list, function(x) x$concordances))

pdf(paste0(figdir, "/summary_orig_vs_mock", exts, dtpext, ".pdf"), 
    width = 10, height = 7)

nbr_keep <- unique(intersect(subset(concs, tp == "original")$ncells1,
                             subset(concs, tp == "mock")$ncells1))
## Remove extension from method name
concs$method <- gsub(exts, "", concs$method)

concsum <- concs %>% as.data.frame() %>%
  dplyr::filter(ncells1 %in% nbr_keep & ncells2 %in% nbr_keep) %>%
  dplyr::filter(ncells1 == ncells2) %>%
  dplyr::group_by(dataset, ncells1, method) %>% 
  dplyr::summarize(tstat = ttest(tp, auc),
                   mediandiff = median(auc[tp == "original"]) - median(auc[tp == "mock"]))

print(concs %>% dplyr::filter(tp == "original") %>%
        dplyr::filter(ncells1 %in% nbr_keep & ncells2 %in% nbr_keep) %>%
        dplyr::filter(ncells1 == ncells2) %>%
        ggplot(aes(x = method, y = auc, col = method)) + 
        geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) + 
        theme_bw() + xlab("") + ylab("area under concordance curve, signal data set") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

print(concs %>% dplyr::filter(tp == "mock") %>%
        dplyr::filter(ncells1 %in% nbr_keep & ncells2 %in% nbr_keep) %>%
        dplyr::filter(ncells1 == ncells2) %>%
        ggplot(aes(x = method, y = auc, col = method)) + 
        geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) + 
        theme_bw() + xlab("") + ylab("area under concordance curve, mock data set") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

print(concs %>% dplyr::filter(tp == "original") %>%
        dplyr::filter(ncells1 %in% nbr_keep & ncells2 %in% nbr_keep) %>%
        dplyr::filter(ncells1 == ncells2) %>%
        ggplot(aes(x = method, y = auc, col = method)) + 
        geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
        geom_point(position = position_jitter(width = 0.2)) + 
        facet_wrap(~dataset) + 
        theme_bw() + xlab("") + ylab("area under concordance curve, signal data set") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

print(concs %>% dplyr::filter(tp == "mock") %>%
        dplyr::filter(ncells1 %in% nbr_keep & ncells2 %in% nbr_keep) %>%
        dplyr::filter(ncells1 == ncells2) %>%
        ggplot(aes(x = method, y = auc, col = method)) + 
        geom_boxplot(outlier.size = -1) + ylim(0, 1) + 
        geom_point(position = position_jitter(width = 0.2)) + 
        facet_wrap(~dataset) + 
        theme_bw() + xlab("") + ylab("area under concordance curve, signal data set") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

print(concsum %>% 
        ggplot(aes(x = method, y = tstat, col = method)) + 
        geom_boxplot(outlier.size = -1) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + ylab("t-statistic, area under concordance curve (original - mock)") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

print(concsum %>% 
        ggplot(aes(x = method, y = mediandiff, col = method)) + 
        geom_boxplot(outlier.size = -1) +
        geom_point(position = position_jitter(width = 0.2), aes(shape = dataset)) +
        theme_bw() + xlab("") + 
        ylab("difference between median area under concordance curve (original - mock)") + 
        scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
        scale_shape_discrete(name = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)))

dev.off()

saveRDS(NULL, paste0(figdir, "/summary_orig_vs_mock", exts, dtpext, ".rds"))

sessionInfo()


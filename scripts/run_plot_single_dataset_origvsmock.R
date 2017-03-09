args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

demethods <- strsplit(demethods, ",")[[1]]

print(demethods)
print(dataset)
print(filt)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")
source("/home/Shared/data/seq/conquer/comparison/scripts/plot_single_dataset_origvsmock.R")

cobras <- list()

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

for (tp in c("", "mock")) {
  cobras[[paste0("tp_", tp)]] <- readRDS(paste0("figures/cobra_data/", dataset, tp, exts, "_cobra.rds"))
  # 
  # ## Create iCOBRA object from the result files for the different methods
  # (resfiles <- paste0("/home/Shared/data/seq/conquer/comparison/results/",
  #                     dataset, tp, "_", demethods, exts, ".rds"))
  # file.exists(resfiles)
  # cobra <- NULL
  # timings <- list()
  # for (rf in resfiles) {
  #   print(rf)
  #   rf <- readRDS(rf)
  #   for (nm in names(rf)) {
  #     print(paste0(nm, ":", paste(names(rf[[nm]]), collapse = ", ")))
  #     timings[[nm]] <- rf[[nm]]$timing
  #     df <- rf[[nm]]$df
  #     if ("pval" %in% colnames(df)) {
  #       cobra <- COBRAData(pval = setNames(data.frame(mt = df$pval,
  #                                                     row.names = rownames(df)), nm),
  #                          object_to_extend = cobra)
  #     }
  #     if ("padj" %in% colnames(df)) {
  #       cobra <- COBRAData(padj = setNames(data.frame(mt = df$padj,
  #                                                     row.names = rownames(df)), nm),
  #                          object_to_extend = cobra)
  #     }
  #     if ("score" %in% colnames(df)) {
  #       cobra <- COBRAData(score = setNames(data.frame(mt = df$score,
  #                                                      row.names = rownames(df)), nm),
  #                          object_to_extend = cobra)
  #     }
  #   }
  # }
  # 
  # cobra <- calculate_adjp(cobra)
  # 
  # cobras[[paste0("tp_", tp)]] <- cobra
}

pdf(paste0("figures/orig_vs_mock/", dataset, exts, "_orig_vs_mock.pdf"), width = 14, height = 9)
summary_data <- plot_compare_orig_mock(cobras, colvec = cols, summary_data = list())
dev.off()

summary_data <- lapply(summary_data, function(L) {
  L$dataset <- dataset
  L$filt <- filt
  L
})

# summary_data$spearman$dataset <- dataset
# summary_data$spearman$filt <- filt
# summary_data$jaccard_adjp0.05$dataset <- dataset
# summary_data$jaccard_adjp0.05$filt <- filt
# summary_data$jaccard_adjp0.1$dataset <- dataset
# summary_data$jaccard_adjp0.1$filt <- filt
# summary_data$jaccard_adjp0.2$dataset <- dataset
# summary_data$jaccard_adjp0.2$filt <- filt

saveRDS(summary_data, file = paste0("figures/orig_vs_mock/", dataset, exts,
                                    "_orig_vs_mock_summary_data.rds"))

sessionInfo()
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(filt)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")
source("/home/Shared/data/seq/conquer/comparison/scripts/plot_single_dataset_origvsmock.R")

concordances <- list()

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

for (tp in c("", "mock")) {
  concordances[[paste0("tp_", tp)]] <- readRDS(paste0("figures/consistency/", 
                                                      dataset, tp, exts, "_concordances.rds"))
}

pdf(paste0("figures/orig_vs_mock/", dataset, exts, "_orig_vs_mock.pdf"), width = 14, height = 9)
summary_data <- plot_compare_orig_mock(concordances, colvec = cols, summary_data = list())
dev.off()

summary_data <- lapply(summary_data, function(L) {
  L$dataset <- dataset
  L$filt <- filt
  L
})

saveRDS(summary_data, file = paste0("figures/orig_vs_mock/", dataset, exts,
                                    "_orig_vs_mock_summary_data.rds"))

sessionInfo()
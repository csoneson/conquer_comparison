args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(filt)
print(concordancedir)
print(figdir)

## Define the number of top-ranked genes to consider for AUC calculations
K0 <- c(100, 1000)

source("scripts/plot_setup.R")
source("scripts/plot_single_dataset_origvsmock.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

concordances <- list()
for (tp in c("", "mock")) {
  concordances[[paste0("tp_", tp)]] <- readRDS(paste0(concordancedir, "/", 
                                                      dataset, tp, exts, "_concordances.rds"))
}

pdf(paste0(figdir, "/", dataset, exts, "_orig_vs_mock.pdf"), width = 14, height = 9)
summary_data <- plot_compare_orig_mock(concordances, colvec = cols, K0 = K0, 
                                       summary_data = list())
dev.off()

summary_data <- lapply(summary_data, function(L) {
  L$dataset <- dataset
  L$filt <- filt
  L
})


saveRDS(summary_data, file = paste0(figdir, "/", dataset, exts,
                                    "_orig_vs_mock_summary_data.rds"))

sessionInfo()
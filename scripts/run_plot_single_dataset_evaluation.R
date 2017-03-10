args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(filt)
print(plottype)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")
source(paste0("/home/Shared/data/seq/conquer/comparison/scripts/plot_single_dataset_", plottype, ".R"))

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)


pdf(paste0("figures/", plottype, "/", dataset, exts, "_", plottype, ".pdf"), width = 16, height = 9)
summary_data <- list()
if (plottype == "timing") {
  timings <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_timings.rds"))
  summary_data <- get(paste0("plot_", plottype))(timings, colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "consistency") {
  cobra <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_cobra.rds"))
  concordances <- readRDS(paste0("figures/consistency/", dataset, exts, "_concordances.rds"))
  summary_data <- get(paste0("plot_", plottype))(cobra, concordances = concordances,
                                                 colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "results_relativetruth_all") {
  cobra <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_cobra.rds"))
  relperf_alltruths <- readRDS(paste0("figures/results_relativetruth_all/", dataset, exts, "_relative_performance.rds"))
  summary_data <- get(paste0("plot_", plottype))(cobra, relperf_alltruths = relperf_alltruths, 
                                                 colvec = cols, exts = exts, summary_data = summary_data)
} else {
  cobra <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_cobra.rds"))
  summary_data <- get(paste0("plot_", plottype))(cobra, colvec = cols, exts = exts, summary_data = summary_data)
}
dev.off()

summary_data <- lapply(summary_data, function(L) {
  L$dataset <- dataset
  L$filt <- filt
  L
})

saveRDS(summary_data, file = paste0("figures/", plottype, "/", dataset, 
                                    exts, "_", plottype, "_summary_data.rds"))

sessionInfo()

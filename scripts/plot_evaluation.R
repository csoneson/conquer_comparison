args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(filt)
print(plottype)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")
source(paste0("/home/Shared/data/seq/conquer/comparison/scripts/plot_", plottype, ".R"))

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

cobra <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_cobra.rds"))
timings <- readRDS(paste0("figures/cobra_data/", dataset, exts, "_timings.rds"))

pdf(paste0("figures/", plottype, "/", dataset, exts, "_", plottype, ".pdf"), width = 14, height = 9)
summary_data <- list()
if (plottype != "timing") {
  summary_data <- get(paste0("plot_", plottype))(cobra, colvec = cols, summary_data = summary_data)
} else {
  summary_data <- get(paste0("plot_", plottype))(timings, colvec = cols, summary_data = summary_data)
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

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datasets)
print(filt)
print(dtpext)
print(summarytype)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")
source(paste0("/home/Shared/data/seq/conquer/comparison/scripts/summarize_", summarytype, ".R"))

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

get(paste0("summarize_", summarytype))(figdir = "figures/summary_crossds", 
                                       datasets = datasets, exts = exts, 
                                       dtpext = dtpext, cols = cols)
saveRDS(NULL, file = paste0("figures/summary_crossds/summary_", summarytype, exts, dtpext, ".rds"))
sessionInfo()


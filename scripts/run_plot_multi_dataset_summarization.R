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
print(singledsfigdir)
print(cobradir)
print(dschardir)
print(concordancedir)
print(figdir)

source("scripts/plot_setup.R")
source(paste0("scripts/summarize_", summarytype, ".R"))

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)

get(paste0("summarize_", summarytype))(figdir = figdir, 
                                       datasets = datasets, exts = exts, 
                                       dtpext = dtpext, cols = cols,
                                       singledsfigdir = singledsfigdir,
                                       cobradir = cobradir,
                                       concordancedir = concordancedir,
                                       dschardir = dschardir)
saveRDS(NULL, file = paste0(figdir, "/summary_", summarytype, exts, dtpext, ".rds"))
sessionInfo()


args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
filt <- strsplit(filt, ",")[[1]]
plotmethods <- strsplit(plotmethods, ",")[[1]]

print(datasets)
print(filt)
print(dtpext)
print(summarytype)
print(plotmethods)
print(singledsfigdir)
print(cobradir)
print(dschardir)
print(concordancedir)
print(origvsmockdir)
print(distrdir)
print(dstypetxt)
print(figdir)

suppressPackageStartupMessages(library(cowplot))
source("scripts/plot_setup.R")
source(paste0("scripts/summarize_", summarytype, ".R"))

if (all(filt == "")) { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
exts <- unique(c("", exts))
cols <- structure(rep(cols, length(exts)), 
                  names = paste0(names(cols), rep(exts, each = length(cols))))

dstypes <- read.delim(dstypetxt, header = TRUE, as.is = TRUE)

plots <- get(paste0("summarize_", summarytype))(figdir = figdir, 
                                                datasets = datasets, exts = exts, 
                                                dtpext = dtpext, cols = cols,
                                                singledsfigdir = singledsfigdir,
                                                cobradir = cobradir,
                                                concordancedir = concordancedir,
                                                dschardir = dschardir,
                                                origvsmockdir = origvsmockdir,
                                                distrdir = distrdir,
                                                plotmethods = plotmethods, 
                                                dstypes = dstypes,
                                                pch_ncells = pch_ncells)
saveRDS(plots, file = paste0(figdir, "/summary_", summarytype, 
                             ifelse(summarytype == "filtering", setdiff(exts, ""), ""), 
                             dtpext, ".rds"))
sessionInfo()


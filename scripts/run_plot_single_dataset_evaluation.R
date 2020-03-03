args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(filt)
print(plottype)
print(cobradir)
print(concordancedir)
print(relperfdir)
print(realperfdir)
print(figdir)

source("scripts/plot_setup.R")
source(paste0("scripts/plot_single_dataset_", plottype, ".R"))

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(names(cols), exts)


pdf(paste0(figdir, "/", plottype, "/", dataset, exts, "_", plottype, ".pdf"), width = 16, height = 9)
summary_data <- list()
if (plottype == "timing") {
  timings <- readRDS(paste0(cobradir, "/", dataset, exts, "_timings.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    timinglist = timings, colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "consistency") {
  cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
  concordances <- readRDS(paste0(concordancedir, "/", dataset, exts, "_concordances.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    cobra = cobra, concordances = concordances,
    colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "results_relativetruth_all") {
  cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
  relperf_alltruths <- readRDS(paste0(relperfdir, "/", 
                                      dataset, exts, "_relative_performance.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    cobra = cobra, relperf_alltruths = relperf_alltruths, 
    colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "performance_realtruth") {
  cobraperf <- readRDS(paste0(realperfdir, "/", dataset, exts, "_performance.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    cobraperf = cobraperf, colvec = cols, exts = exts, summary_data = summary_data)
} else if (plottype == "runfailure") {
  runstatus <- readRDS(paste0(cobradir, "/", dataset, exts, "_runstatus.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    runstatuslist = runstatus, colvec = cols, exts = exts, summary_data = summary_data)
} else {
  cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
  summary_data <- get(paste0("plot_", plottype))(
    cobra = cobra, colvec = cols, exts = exts, summary_data = summary_data)
}
dev.off()

summary_data <- lapply(summary_data, function(L) {
  L$dataset <- dataset
  L$filt <- filt
  L
})

saveRDS(summary_data, file = paste0(figdir, "/", plottype, "/", dataset, 
                                    exts, "_", plottype, "_summary_data.rds"))

sessionInfo()

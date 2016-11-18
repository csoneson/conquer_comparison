args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(config_file)
print(demethod)

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")
source(paste0("/home/Shared/data/seq/conquer/comparison/scripts/apply_", demethod, ".R"))

config <- fromJSON(file = config_file)

print(config)

mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)

subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition

res <- list()

sizes <- names(keep_samples)
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae, keep_samples, sz, i, imposed_condition)
    pdf(paste0(config$figfilebase, "_", demethod, "_", sz, "_", i, ".pdf"))
    res[[paste0(demethod, ".", sz, ".", i)]] <- get(paste0("run_", demethod))(L)
    dev.off()
  }
}

saveRDS(res, file = paste0(config$resfilebase, "_", demethod, ".rds"))

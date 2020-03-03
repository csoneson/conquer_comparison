args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(config_file)
print(demethod)

## filt should be a string of the form type_value_nbr, where type is either TPM 
## or count, and features will be filtered out unless they show at least a 
## "type" of "value" in at least "nbr" samples (if nbr is a single number) or at
## least nbr % of the samples (if nbr ends with p). For now, both value and nbr
## need to be integers (no periods in the name).
print(filt)
(exts <- ifelse(filt == "", "", paste0("_", filt)))

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
source("scripts/prepare_mae.R")
source(paste0("scripts/apply_", demethod, ".R"))

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
    L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i,
                    imposed_condition = imposed_condition, filt = filt, impute = config$impute)
    message(nrow(L$count), " genes, ", ncol(L$count), " samples.")
    pdf(paste0(config$figfilebase, "_", demethod, exts, "_", sz, "_", i, ".pdf"))
    res[[paste0(demethod, exts, ".", sz, ".", i)]] <- get(paste0("run_", demethod))(L)
    res[[paste0(demethod, exts, ".", sz, ".", i)]][["nimp"]] <- L$nimp
    dev.off()
  }
}

saveRDS(res, file = paste0(config$resfilebase, "_", demethod, exts, ".rds"))

sessionInfo()
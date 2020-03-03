## Extend the MultiAssayExperiment for GSE62270-GPL17021 with slots required for
## the analysis

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))

mae <- readRDS("data/GSE62270-GPL17021-orig.rds")

pData(mae) <- mae@colData

keep_cells <- colnames(experiments(mae)[["gene"]])[colSums(assays(experiments(mae)[["gene"]])[["count"]]) > 1000]
mae <- mae[, keep_cells]

assays(experiments(mae)[["gene"]])[["count_lstpm"]] <- assays(experiments(mae)[["gene"]])[["count"]]
assays(experiments(mae)[["gene"]])[["TPM"]] <- 
  sweep(assays(experiments(mae)[["gene"]])[["count"]], 2, 
        colSums(assays(experiments(mae)[["gene"]])[["count"]])/1e6, "/")

summary(colSums(assays(experiments(mae)[["gene"]])[["count_lstpm"]]))
summary(colSums(assays(experiments(mae)[["gene"]])[["TPM"]]))
table(pData(mae)$source_name_ch1)

## Remove transcript counts
experiments(mae)[["tx"]] <- NULL

saveRDS(mae, file = "data/GSE62270-GPL17021.rds")
saveRDS(mae, file = "data/GSE62270-GPL17021mock.rds")

sessionInfo()
date()

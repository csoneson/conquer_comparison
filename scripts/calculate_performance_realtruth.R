args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))

print(dataset)  ## Data set
print(filt)  ## Filtering
print(cobradir)  ## Directory where to look for cobra object (output from prepare_cobra_for_evaluation.R)
print(outdir) ## Directory where to save the calculated concordances

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
## Set p-values and adjusted p-values for untested genes and genes with NA
## values to 1
pval(cobra)[is.na(pval(cobra))] <- 1
padj(cobra)[is.na(padj(cobra))] <- 1
## Similarly, set score for untested genes and genes with NA values to a value
## below the smallest observed score
if (nrow(iCOBRA::score(cobra)) > 0) {
  iCOBRA::score(cobra)[is.na(iCOBRA::score(cobra))] <- min(iCOBRA::score(cobra), na.rm = TRUE) - 1
}

truth <- readRDS(paste0("data/", dataset, "_truth.rds"))

cobra <- COBRAData(truth = truth, 
                   object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "status", 
                                   aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr", "roc"), 
                                   thrs = c(0.01, 0.05, 0.1))

roc(cobraperf) <- roc(cobraperf) %>% dplyr::group_by(method) %>% 
  dplyr::mutate(FPR = c(0, FPR[-1])) %>%
  dplyr::mutate(dFPR = c(0, diff(FPR)),
                dTPR = c(0, diff(TPR)),
                TPRs = c(0, TPR[-length(TPR)])) %>%
  dplyr::mutate(AUC = cumsum(dFPR * dTPR/2 + dFPR * TPRs)) %>%
  dplyr::ungroup() %>% as.data.frame()

saveRDS(cobraperf, file = paste0(outdir, "/", dataset, exts, "_performance.rds"))

sessionInfo()
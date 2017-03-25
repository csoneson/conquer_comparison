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
print(outdir) ## Directory where to save the calculated performances

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

## --------------------------- Data preparation ----------------------------- ##
cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
## Set p-values and adjusted p-values for untested genes and genes with NA
## values to 1
pval(cobra)[is.na(pval(cobra))] <- 1
padj(cobra)[is.na(padj(cobra))] <- 1
## Similarly, set score for untested genes and genes with NA values to a value
## below the smallest observed score
if (nrow(iCOBRA::score(cobra)) > 0) {
  iCOBRA::score(cobra)[is.na(iCOBRA::score(cobra))] <- 
    min(iCOBRA::score(cobra), na.rm = TRUE) - 1
}

truth <- readRDS(paste0("data/", dataset, "_truth.rds"))

cobra <- COBRAData(truth = truth, 
                   object_to_extend = cobra)

## --------------------- Performance calculation ---------------------------- ##
## Calculate performance separately for each data set instance, since the truth 
## is different. Use only genes that are actually tested as the basis for 
## evaluations (exclude the ones that are filtered out in advance).
(szi <- unique(paste0(get_nsamples(colnames(padj(cobra))), ".", 
                      get_repl(colnames(padj(cobra))))))
cobraperf <- NULL
for (s in szi) {
  message(s)
  ## Subset cobra object to given data set instance
  cobratmp <- COBRAData(
    padj = padj(cobra)[, grep(s, colnames(padj(cobra)))],
    score = iCOBRA::score(cobra)[, grep(s, colnames(iCOBRA::score(cobra)))]
    truth = truth(cobra)[, c("gene", "status", 
                             grep(s, colnames(truth(cobra)), value = TRUE))])
  
  ## Keep only genes that are tested for this data set instance
  truth(cobratmp) <- 
    truth(cobratmp)[!is.na(truth(cobratmp)[, paste0("tested.", s)]), , drop = FALSE]
  padj(cobratmp) <- 
    padj(cobratmp)[rownames(padj(cobratmp)) %in% rownames(truth(cobratmp)), , drop = FALSE]
  iCOBRA::score(cobratmp) <- 
    iCOBRA::score(cobratmp)[rownames(iCOBRA::score(cobratmp)) %in% 
                              rownames(truth(cobratmp)), , drop = FALSE]
  
  ## Calculate performance
  cobraperftmp <- calculate_performance(cobratmp, binary_truth = "status", 
                                        aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr", "roc"), 
                                        thrs = c(0.01, 0.05, 0.1))
  if (is.null(cobraperf)) {
    cobraperf <- cobraperftmp
  } else {
    fdrtpr(cobraperf) <- rbind(fdrtpr(cobraperf), fdrtpr(cobraperftmp))
    fdrtprcurve(cobraperf) <- rbind(fdrtprcurve(cobraperf), fdrtprcurve(cobraperftmp))
    fdrnbrcurve(cobraperf) <- rbind(fdrnbrcurve(cobraperf), fdrnbrcurve(cobraperftmp))
    tpr(cobraperf) <- rbind(tpr(cobraperf), tpr(cobraperftmp))
    fpr(cobraperf) <- rbind(fpr(cobraperf), fpr(cobraperftmp))
    roc(cobraperf) <- rbind(roc(cobraperf), roc(cobraperftmp))
  }
}

## Calculate partial (cumulative) AUROCs
roc(cobraperf) <- roc(cobraperf) %>% dplyr::group_by(method) %>% 
  dplyr::mutate(FPR = c(0, FPR[-1])) %>%
  dplyr::mutate(dFPR = c(0, diff(FPR)),
                dTPR = c(0, diff(TPR)),
                TPRs = c(0, TPR[-length(TPR)])) %>%
  dplyr::mutate(AUC = cumsum(dFPR * dTPR/2 + dFPR * TPRs)) %>%
  dplyr::ungroup() %>% as.data.frame()

saveRDS(cobraperf, file = paste0(outdir, "/", dataset, exts, "_performance.rds"))

sessionInfo()
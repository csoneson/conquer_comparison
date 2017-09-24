args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(IHW))

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

## If a method has both score and p-value, remove the score to ascertain that
## ROC curves will (if possible) be calculated from p-values
iCOBRA::score(cobra)[, intersect(colnames(iCOBRA::score(cobra)),
                                 colnames(pval(cobra)))] <- NULL
if (ncol(iCOBRA::score(cobra)) == 0) {
  iCOBRA::score(cobra) <- data.frame()
}

truth <- readRDS(paste0("data/", dataset, "_truth.rds"))

cobra <- COBRAData(truth = truth, 
                   object_to_extend = cobra)

## ----------------- Alternative p-value adjustment ------------------------- ##
## IHW
cobra_ihw <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
pvals <- pval(cobra_ihw)
avecpm <- truth(cobra_ihw)[, grep("avecpm", colnames(truth(cobra_ihw))), drop = FALSE]
cn <- structure(colnames(pvals), names = colnames(pvals))
padjs <- Reduce(function(...) dplyr::full_join(..., by = "gene"), lapply(cn, function(ci) {
  acidx <- which(paste0(get_nsamples(colnames(avecpm)), ".", get_repl(colnames(avecpm))) == 
                   paste0(get_nsamples(ci), ".", get_repl(ci)))
  tmp <- dplyr::full_join(data.frame(pvals[, ci, drop = FALSE]) %>% tibble::rownames_to_column("gene"),
                          data.frame(avecpm[, acidx, drop = FALSE]) %>% tibble::rownames_to_column("gene"))
  tmp <- tmp[rowSums(is.na(tmp)) == 0, , drop = FALSE]
  ihw <- IHW::ihw(pvalues = tmp[, ci], covariates = tmp[, grep("avecpm", colnames(tmp))], alpha = 0.05)
  data.frame(gene = tmp$gene, ihw = IHW::adj_pvalues(ihw), stringsAsFactors = FALSE) %>% setNames(c("gene", ci))
}))
rownames(padjs) <- padjs$gene
padjs$gene <- NULL
padj(cobra_ihw) <- padjs

## Set p-values and adjusted p-values for untested genes and genes with NA
## values to 1
pval(cobra_ihw)[is.na(pval(cobra_ihw))] <- 1
padj(cobra_ihw)[is.na(padj(cobra_ihw))] <- 1

truth <- readRDS(paste0("data/", dataset, "_truth.rds"))

cobra_ihw <- COBRAData(truth = truth, 
                       object_to_extend = cobra_ihw)

## --------------------- Performance calculation ---------------------------- ##
## Calculate performance separately for each data set instance, since the truth 
## is different. Use only genes that are actually tested as the basis for 
## evaluations (exclude the ones that are filtered out in advance).
(szi <- unique(paste0(get_nsamples(colnames(padj(cobra))), ".", 
                      get_repl(colnames(padj(cobra))))))
cobraperf <- NULL
cobraperf_ihw <- NULL
for (s in szi) {
  message(s)
  ## Subset cobra objects to given data set instance
  cobratmp <- COBRAData(
    pval = pval(cobra)[, grep(s, colnames(pval(cobra))), drop = FALSE],
    padj = padj(cobra)[, grep(s, colnames(padj(cobra))), drop = FALSE],
    score = iCOBRA::score(cobra)[, grep(s, colnames(iCOBRA::score(cobra))), drop = FALSE],
    truth = truth(cobra)[, c("gene", "status", 
                             grep(s, colnames(truth(cobra)), value = TRUE)), drop = FALSE])
  
  cobratmp_ihw <- COBRAData(
    pval = pval(cobra_ihw)[, grep(s, colnames(pval(cobra_ihw))), drop = FALSE],
    padj = padj(cobra_ihw)[, grep(s, colnames(padj(cobra_ihw))), drop = FALSE],
    score = iCOBRA::score(cobra_ihw)[, grep(s, colnames(iCOBRA::score(cobra_ihw))), drop = FALSE],
    truth = truth(cobra_ihw)[, c("gene", "status", 
                                 grep(s, colnames(truth(cobra_ihw)), value = TRUE)), drop = FALSE])
  
  ## Keep only genes that are tested for this data set instance
  truth(cobratmp) <- 
    truth(cobratmp)[!is.na(truth(cobratmp)[, paste0("tested.", s)]), , drop = FALSE]
  pval(cobratmp) <- 
    pval(cobratmp)[rownames(pval(cobratmp)) %in% rownames(truth(cobratmp)), , drop = FALSE]
  padj(cobratmp) <- 
    padj(cobratmp)[rownames(padj(cobratmp)) %in% rownames(truth(cobratmp)), , drop = FALSE]
  iCOBRA::score(cobratmp) <- 
    iCOBRA::score(cobratmp)[rownames(iCOBRA::score(cobratmp)) %in% 
                              rownames(truth(cobratmp)), , drop = FALSE]
  
  truth(cobratmp_ihw) <- 
    truth(cobratmp_ihw)[!is.na(truth(cobratmp_ihw)[, paste0("tested.", s)]), , drop = FALSE]
  pval(cobratmp_ihw) <- 
    pval(cobratmp_ihw)[rownames(pval(cobratmp_ihw)) %in% rownames(truth(cobratmp_ihw)), , drop = FALSE]
  padj(cobratmp_ihw) <- 
    padj(cobratmp_ihw)[rownames(padj(cobratmp_ihw)) %in% rownames(truth(cobratmp_ihw)), , drop = FALSE]
  iCOBRA::score(cobratmp_ihw) <- 
    iCOBRA::score(cobratmp_ihw)[rownames(iCOBRA::score(cobratmp_ihw)) %in% 
                                  rownames(truth(cobratmp_ihw)), , drop = FALSE]
  
  ## Calculate performance
  cobraperftmp <- calculate_performance(cobratmp, binary_truth = "status", 
                                        aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr", "roc"), 
                                        thrs = c(0.01, 0.05, 0.1))
  
  cobraperftmp_ihw <- calculate_performance(cobratmp_ihw, binary_truth = "status", 
                                            aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr"), 
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
  
  if (is.null(cobraperf_ihw)) {
    cobraperf_ihw <- cobraperftmp_ihw
  } else {
    fdrtpr(cobraperf_ihw) <- rbind(fdrtpr(cobraperf_ihw), fdrtpr(cobraperftmp_ihw))
    fdrtprcurve(cobraperf_ihw) <- rbind(fdrtprcurve(cobraperf_ihw), fdrtprcurve(cobraperftmp_ihw))
    fdrnbrcurve(cobraperf_ihw) <- rbind(fdrnbrcurve(cobraperf_ihw), fdrnbrcurve(cobraperftmp_ihw))
    tpr(cobraperf_ihw) <- rbind(tpr(cobraperf_ihw), tpr(cobraperftmp_ihw))
    fpr(cobraperf_ihw) <- rbind(fpr(cobraperf_ihw), fpr(cobraperftmp_ihw))
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

saveRDS(list(cobraperf = cobraperf, cobraperf_ihw = cobraperf_ihw), 
        file = paste0(outdir, "/", dataset, exts, "_performance.rds"))

sessionInfo()
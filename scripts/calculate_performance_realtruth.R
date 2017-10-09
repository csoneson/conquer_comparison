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
padjs <- Reduce(function(...) dplyr::full_join(..., by = "gene"), 
                lapply(cn, function(ci) {
                  acidx <- which(paste0(get_nsamples(colnames(avecpm)), ".", 
                                        get_repl(colnames(avecpm))) == 
                                   paste0(get_nsamples(ci), ".", get_repl(ci)))
                  tmp <- dplyr::full_join(data.frame(pvals[, ci, drop = FALSE]) %>% 
                                            tibble::rownames_to_column("gene"),
                                          data.frame(avecpm[, acidx, drop = FALSE]) %>% 
                                            tibble::rownames_to_column("gene"))
                  tmp <- tmp[rowSums(is.na(tmp)) == 0, , drop = FALSE]
                  ihw <- IHW::ihw(pvalues = tmp[, ci], 
                                  covariates = tmp[, grep("avecpm", colnames(tmp))], alpha = 0.05)
                  data.frame(gene = tmp$gene, ihw = IHW::adj_pvalues(ihw), 
                             stringsAsFactors = FALSE) %>% setNames(c("gene", ci))
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

cobras <- list(orig = cobra, ihw = cobra_ihw)
cobraperfs <- list(orig = NULL, ihw = NULL)
for (s in szi) {
  message(s)
  for (nm in names(cobras)) {
    cobratmp <- COBRAData(
      pval = pval(cobras[[nm]])[, grep(s, colnames(pval(cobras[[nm]]))), drop = FALSE],
      padj = padj(cobras[[nm]])[, grep(s, colnames(padj(cobras[[nm]]))), drop = FALSE],
      score = iCOBRA::score(cobras[[nm]])[, grep(s, colnames(iCOBRA::score(cobras[[nm]]))), drop = FALSE],
      truth = truth(cobras[[nm]])[, c("gene", "status", 
                                      grep(s, colnames(truth(cobras[[nm]])), value = TRUE)), drop = FALSE]
    )
    
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
    
    ## Calculate performance
    cobraperftmp <- calculate_performance(cobratmp, binary_truth = "status", 
                                          aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr", "roc"), 
                                          thrs = c(0.01, 0.05, 0.1))
    
    if (is.null(cobraperfs[[nm]])) {
      cobraperfs[[nm]] <- cobraperftmp
    } else {
      fdrtpr(cobraperfs[[nm]]) <- rbind(fdrtpr(cobraperfs[[nm]]), fdrtpr(cobraperftmp))
      fdrtprcurve(cobraperfs[[nm]]) <- rbind(fdrtprcurve(cobraperfs[[nm]]), fdrtprcurve(cobraperftmp))
      fdrnbrcurve(cobraperfs[[nm]]) <- rbind(fdrnbrcurve(cobraperfs[[nm]]), fdrnbrcurve(cobraperftmp))
      tpr(cobraperfs[[nm]]) <- rbind(tpr(cobraperfs[[nm]]), tpr(cobraperftmp))
      fpr(cobraperfs[[nm]]) <- rbind(fpr(cobraperfs[[nm]]), fpr(cobraperftmp))
      roc(cobraperfs[[nm]]) <- rbind(roc(cobraperfs[[nm]]), roc(cobraperftmp))
    }
  }
}
cobraperf <- cobraperfs[["orig"]]
cobraperf_ihw <- cobraperfs[["ihw"]]

## Calculate partial (cumulative) AUROCs
roc(cobraperf) <- roc(cobraperf) %>% dplyr::group_by(method) %>% 
  dplyr::mutate(FPR = c(0, FPR[-1])) %>%
  dplyr::mutate(dFPR = c(0, diff(FPR)),
                dTPR = c(0, diff(TPR)),
                TPRs = c(0, TPR[-length(TPR)])) %>%
  dplyr::mutate(AUC = cumsum(dFPR * dTPR/2 + dFPR * TPRs),
                AUCflat = cumsum(dFPR * TPRs)) %>%
  dplyr::ungroup() %>% as.data.frame()

saveRDS(list(cobraperf = cobraperf, cobraperf_ihw = cobraperf_ihw), 
        file = paste0(outdir, "/", dataset, exts, "_performance.rds"))

sessionInfo()
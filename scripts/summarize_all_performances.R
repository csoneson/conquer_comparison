args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))

print(trueperformancerds)
print(truefprrds)
print(timingrds)
print(fprbiasrds)
print(failureraterds)
print(origvsmockrds)
print(outrds)

## ========================================================================== ##
## Failure rate
## ========================================================================== ##
## Evaluate across all data sets.
## Criteria:
## - Good: average failure rate < 0.01
## - Intermediate: 0.01 < failure rate < 0.25
## - Bad: failure rate > 0.25
##
failurerate <- readRDS(failureraterds)
getfailurerate <- function(failurerate) {
  if (failurerate < 0.01) "good"
  else if (failurerate < 0.25) "intermediate"
  else "bad"
}
failurerate <- data.frame(method = rownames(failurerate),
                          FailureRate = sapply(apply(failurerate, 1, mean), getfailurerate),
                          stringsAsFactors = FALSE)

## ========================================================================== ##
## FDR control
## ========================================================================== ##
## Evaluate after filtering.
## Criteria, "average" FDR control:
## - Good: No more than 75% of FDPs on one side of 0.05, 0.0167 < median FDP < 0.15
## - Intermediate: 0.15 < median FDP < 0.25 or 0.01 < median FDP < 0.0167, or
##   0.0167 < median FDP < 0.15 but more than 75% of FDPs on one side of 0.05
## - Bad: median FDP > 0.25 or < 0.01
##
## Criteria, "worst-case" FDR control:
## - Good: |log2(maximal FDR/0.05)| < 2
## - Intermediate: 2 < |log2(maximal FDR/0.05)| < 2.807
## - Bad: |log2(maximal FDR/0.05)| > 2.807
##
getfdrcontrols <- function(medianfdr, fracout) {
  stopifnot(length(medianfdr) == length(fracout))
  sapply(seq_len(length(medianfdr)), function(i) {
    if (medianfdr[i] < 0.15 && medianfdr[i] > 0.0167 && fracout[i] < 0.75) "good"
    else if (medianfdr[i] < 0.25 && medianfdr[i] > 0.01) "intermediate"
    else "bad"
  })
}
getextrfdr <- function(maxfdr) {
  sapply(seq_len(length(maxfdr)), function(i) {
    if (maxfdr[i] < 0.15) "good"
    else if (maxfdr[i] < 0.35) "intermediate"
    else "bad"
  })
}
trueperf <- readRDS(trueperformancerds)
fdrcontrol <- trueperf %>% dplyr::filter(filt == "TPM_1_25p") %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(frachigh = mean(FDR > 0.05, na.rm = TRUE),
                   fraclow = mean(FDR < 0.05, na.rm = TRUE),
                   medianFDR = median(FDR, na.rm = TRUE),
                   minfdr = min(FDR, na.rm = TRUE),
                   maxfdr = max(FDR, na.rm = TRUE)) %>%
  dplyr::mutate(fracout = pmax(frachigh, fraclow)) %>%
  dplyr::mutate(MedianFDP = getfdrcontrols(medianFDR, fracout),
                MaxFDP = getextrfdr(maxfdr)) %>%
  dplyr::select(method, MedianFDP, MaxFDP)

## ========================================================================== ##
## Power
## ========================================================================== ##
## Evaluate after filtering, and only in data set instances with more than 20 cells
## Criteria:
## - Good: median TPR > 0.8
## - Intermediate: 0.6 < median TPR < 0.8
## - Bad: median TPR < 0.6
##
gettpr <- function(tpr) {
  if (tpr > 0.8) "good"
  else if (tpr > 0.6) "intermediate"
  else "bad"
}
trueperf <- readRDS(trueperformancerds)
power <- trueperf %>% dplyr::filter(filt == "TPM_1_25p") %>%
  dplyr::filter(as.numeric(as.character(ncells_fact)) > 20) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(medianTPR = median(TPR)) %>%
  dplyr::mutate(TPR = sapply(medianTPR, gettpr)) %>%
  dplyr::select(method, TPR)

## ========================================================================== ##
## Area under ROC curve
## ========================================================================== ##
## Evaluate after filtering
## Criteria:
## - Good: median AUC > 0.8
## - Intermediate: 0.65 < median AUC < 0.8
## - Bad: median AUC < 0.65
##
getauc <- function(auc) {
  if (auc > 0.8) "good"
  else if (auc > 0.65) "intermediate"
  else "bad"
}
trueperf <- readRDS(trueperformancerds)
auroc <- trueperf %>% dplyr::filter(filt == "TPM_1_25p") %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(medianAUC = median(AUROC)) %>%
  dplyr::mutate(AUROC = sapply(medianAUC, getauc)) %>%
  dplyr::select(method, AUROC)

## ========================================================================== ##
## Type I error
## ========================================================================== ##
## Evaluate after filtering
## Criteria, "average" type I error control:
## - Good: |log2(median FPR/0.05)| < 1
## - Intermediate: 1 < |log2(median FPR/0.05)| < 2
## - Bad: |log2(median FPR/0.05)| > 2
##
## Criteria, "worst-case" type I error control:
## - Good: |log2(maximal FPR/0.05)| < 1
## - Intermediate: 1 < |log2(maximal FPR/0.05)| < 2.322
## - Bad: |log2(maximal FPR/0.05)| > 2.322
##
fpr <- readRDS(truefprrds)
getfpr <- function(fpr) {
  if (abs(log2(fpr/0.05)) < log2(1.5)) "good"
  else if (abs(log2(fpr/0.05)) < 2) "intermediate"
  else "bad"
}
getextrfpr <- function(maxfpr, minfpr) {
  sapply(seq_len(length(maxfpr)), function(i) {
    if (maxfpr[i] < 0.1) "good"
    else if (maxfpr[i] < 0.25) "intermediate"
    else "bad"
  })
}
fpr <- fpr %>% dplyr::filter(filt == "TPM_1_25p") %>%
  dplyr::group_by(method, dtype) %>%
  dplyr::summarise(medianFPR = median(FPR),
                   maxFPR = max(FPR),
                   minFPR = min(FPR)) %>%
  dplyr::mutate(MedianFPR = sapply(medianFPR, getfpr),
                MaxFPR = getextrfpr(maxFPR, minFPR))
fprumi <- fpr %>% dplyr::filter(dtype == "UMI") %>%
  dplyr::select(method, MedianFPR, MaxFPR) %>%
  dplyr::rename(MedianFPR_UMI = MedianFPR,
                MaxFPR_UMI = MaxFPR)
fprfl <- fpr %>% dplyr::filter(dtype == "full-length") %>%
  dplyr::select(method, MedianFPR, MaxFPR) %>%
  dplyr::rename(MedianFPR_FL = MedianFPR,
                MaxFPR_FL = MaxFPR)

## ========================================================================== ##
## Timing
## ========================================================================== ##
## Evaluate based on all data sets
## Criteria, speed:
## - Good: median relative computation time requirement < 0.1
## - Intermediate: 0.1 < median relative computation time requirement < 0.8
## - Bad: median relative computation time requirement > 0.8
##
## Criteria, scaling with number of cells:
## - Good: median exponent < 0.5
## - Intermediate: 0.5 < median exponent < 1
## - Bad: median exponent > 1
##
timing <- readRDS(timingrds)
getexp <- function(exponent) {
  if (exponent < 0.5) "good"
  else if (exponent < 1) "intermediate"
  else "bad"
}
scaling <- timing$timing_ncelldep %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(medianExponent = median(expn)) %>%
  dplyr::mutate(Scalability = sapply(medianExponent, getexp)) %>%
  dplyr::select(method, Scalability)

getspeed <- function(reltime) {
  if (reltime < 0.1) "good"
  else if (reltime < 0.7) "intermediate"
  else "bad"
}
speed <- timing$timing %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(medianSpeed = median(rel_timing)) %>%
  dplyr::mutate(Speed = sapply(medianSpeed, getspeed)) %>%
  dplyr::select(method, Speed)

## ========================================================================== ##
## FP biases
## ========================================================================== ##
## Evaluate on unfiltered null data
## Criteria:
## - Good: No false positive genes detected, or |median SNR| < 0.5 for all four
##   statistics
## - Intermediate: For at least one statistic, |median SNR| > 0.5, but for all
##   statistics, |median SNR| < 1
## - Bad: For at least one statistics, |median SNR| > 1
##
biases <- readRDS(fprbiasrds)
getsnr <- function(snrs) {
  if (all(abs(snrs[!is.na(snrs)]) < 0.5) | all(is.na(snrs))) "good"
  else if (all(abs(snrs[!is.na(snrs)]) < 1)) "intermediate"
  else "bad"
}
bias <- biases %>%
  dplyr::group_by(method, charac) %>%
  dplyr::summarise(mediansnr = median(snr, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(BiasDEG = getsnr(mediansnr))

## ========================================================================== ##
## Robustness
## ========================================================================== ##
## Evaluate after filtering
## Criteria:
## - Good: The t-statistic of robustness values between signal and null data
##   sets is > 2 for GSE60749 and 10XMonoCytoT, and all t-statistics are > 0
## - Intermediate: Any of the t-statistics for GSE60749 or 10XMonoCytoT is < 2, but
##   all t-statistics are > 0
## - Bad: Any t-statistics is < 0
##
ttest <- function(x, y) {
  tryCatch({
    t.test(x, y, var.equal = FALSE)$stat
  }, error = function(e) NA)
}
consistency <- readRDS(origvsmockrds)
getconsist1 <- function(consist) {
  if (all(consist[!is.na(consist)] > 2)) "good"
  else "intermediate"
}
getconsist2 <- function(consist) {
  if (any(consist[!is.na(consist)] < 0)) "bad"
  else "good"
}
consistency <- consistency %>%
  dplyr::filter(filt == "TPM_1_25p") %>%
  dplyr::filter(k == 100) %>%
  dplyr::group_by(method, dataset, filt, ncells, k) %>%
  dplyr::mutate(tokeep = length(unique(tp)) == 2) %>%
  dplyr::filter(tokeep) %>% 
  dplyr::summarize(tst = ttest(AUCs[tp == "signal"], AUCs[tp == "mock"]))

tmp <- consistency %>%
  dplyr::filter(dataset %in% c("10XMonoCytoT", "GSE60749-GPL13112")) %>%
  dplyr::group_by(method) %>%
  dplyr::summarize(consist1 = getconsist1(tst))
tmp2 <- consistency %>%
  dplyr::group_by(method) %>%
  dplyr::summarize(consist2 = getconsist2(tst))
consist <- dplyr::full_join(tmp, tmp2)
consist$consist1[consist$consist2 == "bad"] <- "bad"
consist <- consist %>% dplyr::rename(Consistency = consist1) %>%
  dplyr::select(method, Consistency)

## ========================================================================== ##
## Complex design
## ========================================================================== ##
## Criteria:
## - Good: Allows arbitrary complex (fixed) designs
## - Intermediate: Can accommodate a limited set of designs
## - Bad: Only performs two-group comparisons
##
complex <- data.frame(ComplexDesign = c(voomlimma = "good", Wilcoxon = "bad",
                                        monoclecensus = "good", limmatrend = "good",
                                        D3E = "bad", edgeRQLF = "good",
                                        MASTtpmDetRate = "good", ROTSvoom = "bad",
                                        MASTcpmDetRate = "good", ttest = "bad",
                                        NODES = "bad", edgeRLRTrobust = "good",
                                        edgeRLRT = "good", MASTtpm = "good",
                                        MASTcpm = "good", SAMseq = "intermediate",
                                        DESeq2 = "good", edgeRLRTdeconv = "good",
                                        edgeRLRTcensus = "good", SeuratTobit = "bad",
                                        metagenomeSeq = "good", monocle = "good",
                                        DESeq2nofilt = "good", DEsingle = "bad",
                                        ROTStpm = "bad", SeuratBimodnofilt = "bad",
                                        BPSC = "good", ROTScpm = "bad",
                                        scDD = "bad", DESeq2census = "good",
                                        SCDE = "bad", SeuratBimod = "bad",
                                        SeuratBimodIsExpr2 = "bad",
                                        monoclecount = "good", 
                                        edgeRQLFDetRate = "good",
                                        zinbwaveedgeR = "good", zinbwaveDESeq2 = "good",
                                        DESeq2betapFALSE = "good", DESeq2LRT = "good",
                                        logregLRT = "bad"),
                      stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column(var = "method")


## ========================================================================== ##
## Summarize
## ========================================================================== ##
allperf <- dplyr::full_join(fdrcontrol, power) %>%
  dplyr::full_join(auroc) %>%
  dplyr::full_join(fprumi) %>%
  dplyr::full_join(fprfl) %>%
  dplyr::full_join(scaling) %>%
  dplyr::full_join(speed) %>%
  dplyr::full_join(bias) %>%
  dplyr::full_join(consist) %>%
  dplyr::left_join(complex) %>%
  dplyr::left_join(failurerate) %>%
  dplyr::mutate(BiasDEG = replace(BiasDEG, is.na(BiasDEG), "good"))

## ========================================================================== ##
## Convert to numeric values and plot
## ========================================================================== ##
allperf[allperf == "good"] <- 2
allperf[allperf == "intermediate"] <- 1
allperf[allperf == "bad"] <- 0
allperf <- as.data.frame(allperf)
tmp <- allperf$method
allperf <- allperf %>% select(-method)
allperf <- allperf %>% mutate_if(is.character, as.numeric)
allperf <- as.data.frame(allperf)
rownames(allperf) <- tmp
## Order by average score
allperf <- allperf[order(rowMeans(allperf, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]

write.table(cbind(method = rownames(allperf), allperf), 
            file = "export_results/Figure5.csv", row.names = FALSE, col.names = TRUE, 
            sep = ",", quote = FALSE)

pdf(gsub("rds$", "pdf", outrds), width = 6, height = 8)
pheatmap(allperf, cluster_rows = FALSE, cluster_cols = FALSE,
         color = c("#D55E00", "#F0E442", "#56B4E9"), 
         #color = c("#E8601C", "#F6C141", "#90C987"), 
         breaks = c(-0.5, 0.5, 1.5, 2.5),
         scale = "none", legend_breaks = c(0, 1, 2),
         legend_labels = c("poor", "intermediate", "good"), fontsize = 11,
         gaps_col = seq_len(ncol(allperf)),
         gaps_row = seq_len(nrow(allperf)))
dev.off()

saveRDS(allperf, file = outrds)
sessionInfo()
date()

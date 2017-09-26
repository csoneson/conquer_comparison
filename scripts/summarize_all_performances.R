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
print(outrds)

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
    else if (medianfdr[i] > 0.25 || medianfdr[i] < 0.01) "bad"
  })
}
getextrfdr <- function(maxfdr) {
  sapply(seq_len(length(maxfdr)), function(i) {
    if (abs(log2(maxfdr[i]/0.05)) < log2(4)) "good"
    else if (abs(log2(maxfdr[i]/0.05)) < log2(7)) "intermediate"
    else "bad"
  })
}
trueperf <- readRDS(trueperformancerds)
fdrcontrol <- trueperf$FDR_all_byfdrcontrol_TPM_1_25p$data %>%
  dplyr::group_by(method) %>% 
  dplyr::summarise(frachigh = mean(FDR > 0.05), 
                   fraclow = mean(FDR < 0.05),
                   medianFDR = median(FDR),
                   minfdr = min(FDR),
                   maxfdr = max(FDR)) %>%
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
power <- trueperf$TPR_all_byfdrcontrol_TPM_1_25p$data %>%
  dplyr::filter(as.numeric(as.character(n_samples)) > 20) %>% 
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
auroc <- trueperf$auroc_all_byfdrcontrol_TPM_1_25p$data %>%
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
    if (abs(log2(maxfpr[i]/0.05)) < log2(2)) "good"
    else if (abs(log2(maxfpr[i]/0.05)) < log2(5)) "intermediate"
    else "bad"
  })
}
fpr <- fpr$truefpr_sep_TPM_1_25p$data %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarise(medianFPR = median(FPR),
                   maxFPR = max(FPR),
                   minFPR = min(FPR)) %>%
  dplyr::mutate(MedianFPR = sapply(medianFPR, getfpr),
                MaxFPR = getextrfpr(maxFPR, minFPR)) %>%
  dplyr::select(method, MedianFPR, MaxFPR)

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
scaling <- timing$timing_exponent_ncells$data %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarise(medianExponent = median(expn)) %>%
  dplyr::mutate(Scalability = sapply(medianExponent, getexp)) %>%
  dplyr::select(method, Scalability)

getspeed <- function(reltime) {
  if (reltime < 0.1) "good"
  else if (reltime < 0.7) "intermediate"
  else "bad"
}
speed <- timing$rel_timing_boxplot_comb_log$data %>% 
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
  if (all(snrs < 0.5)) "good"
  else if (all(snrs < 1)) "intermediate"
  else "bad"
}
bias <- biases$snr_bystat_$data %>% 
  dplyr::group_by(method, charac) %>%
  dplyr::summarise(mediansnr = median(snr)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(BiasDEG = getsnr(mediansnr))

## ========================================================================== ##
## Robustness
## ========================================================================== ##
## Evaluate after filtering
## Criteria:
## - Good: The t-statistic of robustness values between signal and null data
##   sets is > 1 for GSE60749 and GSE45719, and all t-statistics are > 0
## - Intermediate: Any of the t-statistics for GSE60749 or GSE45719 is < 1, but
##   all t-statistics are > 0
## - Bad: Any t-statistics is < 0
##
consistency <- readRDS("figures/multi_dataset/orig_vs_mock/summary_orig_vs_mock_real.rds")
getconsist1 <- function(consist) {
  if (all(consist[!is.na(consist)] > 2)) "good"
  else "intermediate"
}
getconsist2 <- function(consist) {
  if (any(consist[!is.na(consist)] < 0)) "bad"
  else "good"
}
tmp <- consistency$tstat_auc_TPM_1_25p_100$data %>% 
  dplyr::filter(dataset %in% c("GSE45719", "GSE60749-GPL13112")) %>% 
  dplyr::group_by(method) %>%
  dplyr::summarize(consist1 = getconsist1(tstat))
tmp2 <- consistency$tstat_auc_TPM_1_25p_100$data %>% 
  dplyr::group_by(method) %>%
  dplyr::summarize(consist2 = getconsist2(tstat))
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
                                        SeuratBimodIsExpr2 = "bad"),
                      stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column(var = "method")


## ========================================================================== ##
## Summarize
## ========================================================================== ##
allperf <- dplyr::full_join(fdrcontrol, power) %>%
  dplyr::full_join(auroc) %>%
  dplyr::full_join(fpr) %>%
  dplyr::full_join(scaling) %>%
  dplyr::full_join(speed) %>%
  dplyr::full_join(bias) %>%
  dplyr::full_join(consist) %>%
  dplyr::full_join(complex) %>% 
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
allperf <- allperf[order(rowMeans(allperf, na.rm = TRUE), decreasing = TRUE), ]

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 8)
pheatmap(allperf, cluster_rows = FALSE, cluster_cols = FALSE, 
         color = c("#E8601C", "#F6C141", "#90C987"), breaks = c(-0.5, 0.5, 1.5, 2.5), 
         scale = "none", legend_breaks = c(0, 1, 2), 
         legend_labels = c("poor", "intermediate", "good"), fontsize = 14, 
         gaps_col = seq_len(ncol(allperf)),
         gaps_row = seq_len(nrow(allperf)))
dev.off()

saveRDS(NULL, file = outrds)
sessionInfo()
date()

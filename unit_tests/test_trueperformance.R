## Test true performance calculation

topdir <- "/home/Shared/data/seq/conquer/comparison"

suppressPackageStartupMessages(library(iCOBRA))

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

test_that("true performance is correctly calculated", {
  gw <- getwd()
  setwd(topdir)
  for (ds in c("GSE45719sim123", "GSE74596sim123", "GSE48968-GPL13112sim123")) {
    for (f in c("", "_TPM_1_25p")) {
      cobra <- readRDS(paste0(topdir, "/output/cobra_data/", ds, f, "_cobra.rds"))
      perf <- readRDS(paste0(topdir, "/output/performance_realtruth/", ds, f, "_performance.rds"))
      if ("cobraperf" %in% names(perf)) {
        perf <- perf$cobraperf
      }
      truth <- readRDS(paste0(topdir, "/data/", ds, "_truth.rds"))

      for (mth in colnames(padj(cobra))) {
        trt <- truth(cobra) %>% 
          dplyr::select_("gene", paste0("tested.", get_nsamples(mth), ".", get_repl(mth))) %>%
          dplyr::rename_("tested" = paste0("tested.", get_nsamples(mth), ".", get_repl(mth))) %>%
          dplyr::filter(tested == TRUE) %>% dplyr::inner_join(truth, by = c("gene" = "id")) %>%
          dplyr::inner_join(data.frame(gene = rownames(padj(cobra)), padj = padj(cobra)[, mth],
                                       stringsAsFactors = FALSE), by = "gene") %>%
          dplyr::mutate(padj = replace(padj, is.na(padj), 1))
        fpr_df <- fpr(perf) %>% dplyr::filter(method == mth)
        tpr_df <- tpr(perf) %>% dplyr::filter(method == mth)
        fdrtpr_df <- fdrtpr(perf) %>% dplyr::filter(method == mth)
        
        expect_equal(fpr_df$TOT_CALLED, rep(nrow(trt), nrow(fpr_df)))
        expect_equal(tpr_df$TOT_CALLED, rep(nrow(trt), nrow(tpr_df)))
        expect_equal(fdrtpr_df$TOT_CALLED, rep(nrow(trt), nrow(fdrtpr_df)))
        
        expect_equal(fpr_df$TP, tpr_df$TP)
        expect_equal(fpr_df$FP, tpr_df$FP)
        expect_equal(fpr_df$TN, tpr_df$TN)
        expect_equal(fpr_df$FN, tpr_df$FN)
        
        expect_equal(fpr_df$TP, fdrtpr_df$TP)
        expect_equal(fpr_df$FP, fdrtpr_df$FP)
        expect_equal(fpr_df$TN, fdrtpr_df$TN)
        expect_equal(fpr_df$FN, fdrtpr_df$FN)
        
        expect_equal(fpr_df$NBR, sapply(fpr_df$thr, function(i) {
          length(which(trt$padj <= as.numeric(gsub("^thr", "", i))))
        }))
        expect_equal(fpr_df$TP, sapply(fpr_df$thr, function(i) {
          length(intersect(which(trt$padj <= as.numeric(gsub("^thr", "", i))),
                           which(trt$status == 1)))
        }))
        expect_equal(fpr_df$FP, sapply(fpr_df$thr, function(i) {
          length(intersect(which(trt$padj <= as.numeric(gsub("^thr", "", i))),
                           which(trt$status == 0)))
        }))
        expect_equal(fpr_df$TN, sapply(fpr_df$thr, function(i) {
          length(intersect(which(trt$padj > as.numeric(gsub("^thr", "", i))),
                           which(trt$status == 0)))
        }))
        expect_equal(fpr_df$FN, sapply(fpr_df$thr, function(i) {
          length(intersect(which(trt$padj > as.numeric(gsub("^thr", "", i))),
                           which(trt$status == 1)))
        }))
        
        expect_equal(fpr_df$FPR, fpr_df$FP/fpr_df$NONDIFF)
        expect_equal(tpr_df$TPR, tpr_df$TP/tpr_df$DIFF)
        tmp <- fdrtpr_df$FP/fdrtpr_df$NBR
        tmp[fdrtpr_df$NBR == 0] <- 0
        expect_equal(fdrtpr_df$FDR, tmp)

        ## To test similarity between values in ROC df and those in TPR df, 
        ## ROC curve must be calculated from adjusted p-values (so that CUTOFF
        ## values correspond to p.adj)
        # roc_df <- roc(perf) %>% dplyr::filter(method == mth)
        # roc_df <- roc_df[match(sapply(fpr_df$thr, function(i) {
        #   max(roc_df$ROC_CUTOFF[round(roc_df$ROC_CUTOFF, 15) <= as.numeric(gsub("^thr", "", i))])
        # }), roc_df$ROC_CUTOFF), ]
        # expect_equal(roc_df$FPR, fpr_df$FPR)
        # expect_equal(roc_df$TPR, tpr_df$TPR)
      }
    }        
  }
  setwd(gw)
})
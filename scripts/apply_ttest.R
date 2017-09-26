suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(genefilter))

run_ttest <- function(L) {
  message("t-test")
  session_info <- sessionInfo()
  timing <- system.time({
    tmm <- edgeR::calcNormFactors(L$tpm)
    tpmtmm <- edgeR::cpm(L$tpm, lib.size = tmm * colSums(L$tpm))
    logtpm <- log2(tpmtmm + 1)
    idx <- seq_len(nrow(logtpm))
    names(idx) <- rownames(logtpm)
    ttest_p <- sapply(idx, function(i) {
      t.test(logtpm[i, ] ~ L$condt)$p.value
    })
  })
  
  hist(ttest_p, 50)
  
  list(session_info = session_info,
       timing = timing,
       df = data.frame(pval = ttest_p,
                       row.names = names(ttest_p)))
}
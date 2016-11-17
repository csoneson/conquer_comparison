suppressPackageStartupMessages(library(edgeR))

run_Wilcoxon <- function(L) {
  message("Wilcoxon")
  session_info <- sessionInfo()
  timing <- system.time({
    tmm <- edgeR::calcNormFactors(L$tpm)
    tpmtmm <- edgeR::cpm(L$tpm, lib.size = tmm * colSums(L$tpm))
    idx <- 1:nrow(tpmtmm)
    names(idx) <- rownames(tpmtmm)
    wilcox_p <- sapply(idx, function(i) {
      wilcox.test(tpmtmm[i, ] ~ L$condt)$p.value
    })
  })
  
  hist(wilcox_p, 50)
  
  list(session_info = session_info,
       timing = timing,
       df = data.frame(pval = wilcox_p,
                       row.names = names(wilcox_p)))
}
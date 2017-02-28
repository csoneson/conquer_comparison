suppressPackageStartupMessages(library(MAST))

run_MASTtpm <- function(L) {
  message("MAST, TPM")
  session_info <- sessionInfo()
  timing <- system.time({
    stopifnot(all(names(L$condt) == colnames(L$tpm)))
    grp <- L$condt
    sca <- FromMatrix(exprsArray = log2(L$tpm + 1), 
                      cData = data.frame(wellKey = names(grp), grp = grp))
    zlmdata <- zlm.SingleCellAssay(~grp, sca)
    mast <- lrTest(zlmdata, "grp")
  })
  
  hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
}
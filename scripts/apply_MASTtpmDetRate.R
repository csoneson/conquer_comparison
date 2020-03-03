suppressPackageStartupMessages(library(MAST))

run_MASTtpmDetRate <- function(L) {
  message("MAST, TPM (including detection rate)")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      stopifnot(all(names(L$condt) == colnames(L$tpm)))
      grp <- L$condt
      cdr <- scale(colMeans(L$tpm > 0))
      sca <- FromMatrix(exprsArray = log2(L$tpm + 1), 
                        cData = data.frame(wellKey = names(grp), 
                                           grp = grp, cdr = cdr))
      zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
      mast <- lrTest(zlmdata, "grp")
    })
    
    hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
    
    list(session_info = session_info,
         timing = timing,
         res = mast,
         df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
  }, error = function(e) {
    "MASTtpmDetRate results could not be calculated"
    list(session_info = session_info)
  })
}


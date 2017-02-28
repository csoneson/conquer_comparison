suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))

run_MASTcpm <- function(L) {
  message("MAST, CPM")
  session_info <- sessionInfo()
  timing <- system.time({
    stopifnot(all(names(L$condt) == colnames(L$count)))
    grp <- L$condt
    dge <- DGEList(counts = L$count)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
    sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                      cData = data.frame(wellKey = names(grp), 
                                         grp = grp))
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
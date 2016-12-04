suppressPackageStartupMessages(library(MAST))

run_MASTcountsDetRate <- function(L) {
  message("MAST, counts (including detection rate")
  session_info <- sessionInfo()
  timing <- system.time({
    grp <- L$condt
    countmelt <- reshape2::melt(as.matrix(L$count))
    colnames(countmelt) <- c("gene", "cell", "value")
    sca <- FromFlatDF(countmelt, idvars = "cell", primerid = "gene", 
                      measurement = "value")
    cdr <- colSums(assay(sca) > 0)
    colData(sca)$grp <- grp[match(rownames(colData(sca)), names(grp))]
    colData(sca)$cngeneson <- scale(cdr)
    zlmdata <- zlm.SingleCellAssay(~ cngeneson + grp, sca)
    mast <- waldTest(zlmdata, CoefficientHypothesis(paste0("grp", levels(factor(grp))[2])))
  })
  
  hist(mast[, , "Pr(>Chisq)"][, "hurdle"], 50)
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = data.frame(pval = mast[, , "Pr(>Chisq)"][, "hurdle"],
                       row.names = names(mast[, , "Pr(>Chisq)"][, "hurdle"])))
}
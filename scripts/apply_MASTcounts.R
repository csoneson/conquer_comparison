suppressPackageStartupMessages(library(MAST))

run_MASTcounts <- function(L) {
  message("MAST, counts")
  session_info <- sessionInfo()
  timing <- system.time({
    grp <- L$condt
    countmelt <- reshape2::melt(as.matrix(L$count))
    colnames(countmelt) <- c("gene", "cell", "value")
    sca <- FromFlatDF(countmelt, idvars = "cell", primerid = "gene", 
                      measurement = "value")
    colData(sca)$grp <- grp[match(rownames(colData(sca)), names(grp))]
    zlmdata <- zlm.SingleCellAssay(~grp, sca)
    mast <- waldTest(zlmdata, CoefficientHypothesis(paste0("grp", levels(factor(grp))[2])))
  })
  
  hist(mast[, , "Pr(>Chisq)"][, "hurdle"], 50)
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = data.frame(pval = mast[, , "Pr(>Chisq)"][, "hurdle"],
                       row.names = names(mast[, , "Pr(>Chisq)"][, "hurdle"])))
}
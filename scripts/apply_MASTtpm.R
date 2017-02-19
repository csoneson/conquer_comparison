suppressPackageStartupMessages(library(MAST))

run_MASTtpm <- function(L) {
  message("MAST, TPMs")
  session_info <- sessionInfo()
  timing <- system.time({
    grp <- L$condt
    sca <- FromMatrix(exprsArray = log2(L$tpm + 1), 
                      cData = data.frame(grp = grp, row.names = names(grp)))
    #tpmmelt <- reshape2::melt(as.matrix(L$tpm))
    #colnames(tpmmelt) <- c("gene", "cell", "value")
    #sca <- FromFlatDF(tpmmelt, idvars = "cell", primerid = "gene", 
    #                  measurement = "value")
    #colData(sca)$grp <- grp[match(rownames(colData(sca)), names(grp))]
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
suppressPackageStartupMessages(library(BPSC))
suppressPackageStartupMessages(library(edgeR))

run_BPSC <- function(L) {
  message("BPSC")
  session_info <- sessionInfo()
  timing <- system.time({
    cpms <- edgeR::cpm(L$count, lib.size = colSums(L$count) * edgeR::calcNormFactors(L$count))
    controlIds <- which(L$condt == levels(factor(L$condt))[1])
    design <- model.matrix(~ L$condt)
    coef <- 2
    resbp <- BPglm(data = cpms, controlIds = controlIds, 
                   design = design, coef = coef, estIntPar = FALSE) 
  })
  
  hist(resbp$PVAL, 50)
  
  list(session_info = session_info,
       timing = timing,
       res = resbp,
       df = data.frame(pval = resbp$PVAL,
                       row.names = rownames(L$count)))
}
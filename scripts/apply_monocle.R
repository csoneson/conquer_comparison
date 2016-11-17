suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(edgeR))

run_monocle <- function(L) {
  message("monocle")
  session_info <- sessionInfo()
  timing <- system.time({
    tmm <- edgeR::calcNormFactors(L$tpm)
    tpmtmm <- edgeR::cpm(L$tpm, lib.size = tmm * colSums(L$tpm))
    mon <- newCellDataSet(as.matrix(tpmtmm), 
                          phenoData = new("AnnotatedDataFrame", 
                                          data = data.frame(condition = L$condt, 
                                                            row.names = colnames(tpmtmm))))
    monres <- differentialGeneTest(mon, fullModelFormulaStr = "expression ~ condition")
  })
  
  hist(monres$pval, 50)
  hist(monres$qval, 50)
  
  list(session_info = session_info,
       timing = timing,
       res = monres,
       df = data.frame(pval = monres$pval, 
                       padj = monres$qval,
                       row.names = rownames(monres)))
}
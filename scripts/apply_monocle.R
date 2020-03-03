suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(edgeR))

run_monocle <- function(L) {
  message("monocle")
  session_info <- sessionInfo()
  timing <- system.time({
    mon <- newCellDataSet(as.matrix(L$tpm), 
                          phenoData = new("AnnotatedDataFrame", 
                                          data = data.frame(condition = L$condt, 
                                                            row.names = colnames(L$tpm))),
                          expressionFamily = tobit())
    monres <- differentialGeneTest(mon, fullModelFormulaStr = " ~ condition")
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
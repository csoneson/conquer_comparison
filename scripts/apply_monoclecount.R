suppressPackageStartupMessages(library(monocle))

run_monoclecount <- function(L) {
  message("monoclecount")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      mon <- newCellDataSet(as.matrix(L$count), 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$count))),
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())
      mon <- estimateSizeFactors(mon)
      mon <- estimateDispersions(mon)
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
  }, error = function(e) {
    "monoclecount results could not be calculated"
    list(session_info = session_info)
  })
}

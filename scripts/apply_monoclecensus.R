suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_monoclecensus <- function(L) {
  message("monoclecensus")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      mon <- newCellDataSet(as.matrix(L$tpm), 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$tpm))),
                            expressionFamily = tobit())
      rpc_matrix <- relative2abs(mon)
      mon <- newCellDataSet(cellData = as.matrix(rpc_matrix), 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$tpm))),
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
    "monoclecensus results could not be calculated"
    list(session_info = session_info)
  })
}

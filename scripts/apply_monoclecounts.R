suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_monoclecounts <- function(L) {
  message("monoclecounts")
  session_info <- sessionInfo()
  timing <- system.time({
    mon <- newCellDataSet(as.matrix(L$count), 
                          phenoData = new("AnnotatedDataFrame", 
                                          data = data.frame(condition = L$condt, 
                                                            row.names = colnames(L$count))),
                          expressionFamily = negbinomial.size(),
                          lowerDetectionLimit = 1)
    mon <- DESeq2::estimateSizeFactors(mon)
    mon <- estimateDispersions(mon, fitType = "local")
    monres <- differentialGeneTest(mon, fullModelFormulaStr = "expression ~ condition",
                                   relative_expr = FALSE)
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

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_DESeq2census <- function(L) {
  message("DESeq2census")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      cds <- newCellDataSet(L$tpm, 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$tpm))))
      censuscounts <- relative2abs(cds)
      dds <- DESeqDataSetFromMatrix(countData = round(censuscounts), 
                                    colData = data.frame(condition = L$condt), 
                                    design = ~condition)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("condition", levels(factor(L$condt))[1], 
                                       levels(factor(L$condt))[2]), alpha = 0.05)
    })
    
    plotDispEsts(dds)
    plotMA(res)
    summary(res)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$pval,
                         padj = res$padj,
                         row.names = rownames(res)))
  }, error = function(e) {
    "DESeq2 results could not be calculated"
    list(session_info = session_info)
  })
}
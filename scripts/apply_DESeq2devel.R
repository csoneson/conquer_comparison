suppressPackageStartupMessages(library(devtools))
devtools::load_all("software/DESeq2")

run_DESeq2devel <- function(L) {
  message("DESeq2devel")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      condition <- L$condt
      design <- model.matrix(~ condition)
      dds <- DESeqDataSetFromMatrix(countData = round(L$count), 
                                    colData = data.frame(condition = L$condt), 
                                    design = design)
      dds <- DESeq(dds)
      res <- results(dds, name = resultsNames(dds)[2], alpha = 0.05)
    })
    
    plotDispEsts(dds)
    plotMA(res)
    summary(res)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$pvalue,
                         padj = res$padj,
                         row.names = rownames(res)))
  }, error = function(e) {
    "DESeq2devel results could not be calculated"
    list(session_info = session_info)
  })
}
suppressPackageStartupMessages(library(DEsingle))

run_DEsingle <- function(L) {
  message("DEsingle")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      stopifnot(all(colnames(L$count) == names(L$condt)))
      results <- DEsingle(counts = round(L$count), group = factor(L$condt))
      results <- data.frame(results, stringsAsFactors = FALSE)
    })
    
    hist(results$pvalue, 50)
    hist(results$pvalue.adj.FDR, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = results,
         df = data.frame(pval = results$pvalue,
                         padj = results$pvalue.adj.FDR,
                         row.names = rownames(results)))
  }, error = function(e) {
    "DEsingle results could not be calculated"
    list(session_info = session_info)
  })
}

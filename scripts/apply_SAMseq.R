suppressPackageStartupMessages(library(samr))

run_SAMseq <- function(L) {
  message("SAMseq")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      SAMseq.test <- SAMseq(round(L$count), as.numeric(as.factor(L$condt)),
                            resp.type = "Two class unpaired",
                            geneid = rownames(L$count), genenames = rownames(L$count),
                            nperms = 100, nresamp = 20, fdr.output = 1)
      SAMseq.result.table <- rbind(SAMseq.test$siggenes.table$genes.up,
                                   SAMseq.test$siggenes.table$genes.lo)
      SAMseq.FDR <- rep(NA, nrow(L$count))
      SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(L$count))] <- 
        as.numeric(SAMseq.result.table[, 5])/100
    })
    
    hist(SAMseq.FDR, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = SAMseq.result.table, 
         df = data.frame(padj = SAMseq.FDR, 
                         row.names = rownames(L$count)))
  }, error = function(e) {
    "SAMseq results could not be calculated"
    list(session_info = session_info)
  })
}
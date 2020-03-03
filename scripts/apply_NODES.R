suppressPackageStartupMessages(library(NODES))

run_NODES <- function(L) {
  message("NODES")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      norm_counts <- pQ(L$count)
      grp <- make.names(as.character(L$condt[match(colnames(norm_counts), names(L$condt))]))
      nodes <- NODES(norm_counts, grp)
    })
    
    hist(nodes$qvalues, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = nodes,
         df = data.frame(padj = nodes$qvalues,
                         row.names = rownames(nodes)))
  }, error = function(e) {
    "NODES results could not be calculated"
    list(session_info = session_info)
  })
}
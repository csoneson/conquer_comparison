suppressPackageStartupMessages(library(Seurat))

run_SeuratTobitnofilt <- function(L) {
  message("SeuratTobit")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      tmptpm <- L$tpm
      colnames(tmptpm) <- paste0(L$condt, "__", 1:ncol(L$tpm))
      seur <- new("seurat", raw.data = tmptpm)
      seur <- Setup(seur, project = "scrnaseq", min.cells = 1, min.genes = 0, 
                    do.logNormalize = FALSE, total.expr = 1e4, names.field = 1, 
                    names.delim = "__")
      res <- FindMarkers(seur, ident.1 = levels(factor(L$condt))[1], 
                         ident.2 = levels(factor(L$condt))[2], thresh.use = -Inf, 
                         test.use = "tobit", min.pct = 0, min.cells = 0)
    })
    
    hist(res$p_val, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$p_val,
                         padj = p.adjust(res$p_val, method = "BH"),
                         row.names = rownames(res)))
  }, error = function(e) {
    "Seurat results could not be calculated"
    list(session_info = session_info)
  })
}
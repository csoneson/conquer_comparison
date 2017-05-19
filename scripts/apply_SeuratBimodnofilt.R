suppressPackageStartupMessages(library(Seurat))

run_SeuratBimodnofilt <- function(L) {
  message("SeuratBimod")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      tmpcount <- L$count
      colnames(tmpcount) <- paste0(L$condt, "__", 1:ncol(L$count))
      seur <- new("seurat", raw.data = tmpcount)
      seur <- Setup(seur, project = "scrnaseq", min.cells = 1, min.genes = 0, 
                    do.logNormalize = TRUE, total.expr = 1e4, names.field = 1, 
                    names.delim = "__")
      res <- FindMarkers(seur, ident.1 = levels(factor(L$condt))[1], 
                         ident.2 = levels(factor(L$condt))[2], thresh.use = -Inf, 
                         test.use = "bimod", min.pct = 0, min.cells = 0)
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
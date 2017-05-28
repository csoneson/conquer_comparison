suppressPackageStartupMessages(library(Seurat))

run_SeuratBimodIsExpr2 <- function(L) {
  message("SeuratBimodIsExpr2")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      tmpcount <- L$count
      colnames(tmpcount) <- paste0(L$condt, "__", 1:ncol(L$count))
      seur <- new("seurat", raw.data = tmpcount)
      seur <- Setup(seur, project = "scrnaseq", do.logNormalize = TRUE, 
                    names.field = 1, names.delim = "__", is.expr = 2)
      res <- FindMarkers(seur, ident.1 = levels(factor(L$condt))[1], 
                         ident.2 = levels(factor(L$condt))[2], test.use = "bimod")
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
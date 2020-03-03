suppressPackageStartupMessages(library(Seurat))

run_SeuratTobit <- function(L) {
  message("SeuratTobit")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      tmptpm <- L$tpm
      colnames(tmptpm) <- paste0(L$condt, "__", 1:ncol(L$tpm))
      seur <- new("seurat", raw.data = tmptpm)
      seur <- Setup(seur, project = "scrnaseq", do.logNormalize = FALSE, 
                    names.field = 1, names.delim = "__")
      res <- FindMarkers(seur, ident.1 = levels(factor(L$condt))[1], 
                         ident.2 = levels(factor(L$condt))[2], test.use = "tobit")
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
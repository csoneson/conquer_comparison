suppressPackageStartupMessages(library(metagenomeSeq))

run_metagenomeSeq <- function(L) {
  message("metagenomeSeq")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      obj <- newMRexperiment(L$count, 
                             phenoData = new("AnnotatedDataFrame", 
                                             data = data.frame(condition = L$condt, 
                                                               row.names = colnames(L$count))))
      p <- cumNormStatFast(obj)
      obj <- cumNorm(obj, p = p)
      mod <- model.matrix(~ condition, data = pData(obj))
      res <- fitFeatureModel(obj = obj, mod = mod, coef = 2)
      tbl <- MRtable(obj = res, number = Inf, by = 2, adjustMethod = "BH")
    })
    
    hist(tbl$pvalues, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = tbl$pvalues,
                         padj = tbl$adjPvalues,
                         row.names = rownames(tbl)))
  }, error = function(e) {
    "metagenomeSeq results could not be calculated"
    list(session_info = session_info)
  })
}
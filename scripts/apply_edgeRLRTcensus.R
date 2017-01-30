suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(monocle))

run_edgeRLRTcensus <- function(L) {
  message("edgeRLRTcensus")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      cds <- newCellDataSet(L$tpm, 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$tpm))))
      censuscounts <- relative2abs(cds)
      dge <- DGEList(censuscounts, group = L$condt)
      dge <- calcNormFactors(dge)
      design <- model.matrix(~L$condt)
      dge <- estimateDisp(dge, design = design)
      fit <- glmFit(dge, design = design)
      lrt <- glmLRT(fit)
      tt <- topTags(lrt, n = Inf)
    })
    
    plotBCV(dge)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
    plotSmear(lrt)
    
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  }, error = function(e) {
    "edgeRLRTcensus results could not be calculated"
    list(session_info = session_info)
  })
}

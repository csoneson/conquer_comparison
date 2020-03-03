suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))

run_edgeRLRTdeconv <- function(L) {
  message("edgeRdeconv")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      sce <- newSCESet(countData = L$count)
      sce <- computeSumFactors(sce, sizes = unique(pmin(c(20, 40, 60, 80, 100), ncol(sce)/2)), 
                               positive = TRUE)
      ##isSpike(sce) <- rep(c(FALSE), nrow(sce))
      dge <- convertTo(sce, type = "edgeR")
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
    "edgeRLRTdeconv results could not be calculated"
    list(session_info = session_info)
  })
}
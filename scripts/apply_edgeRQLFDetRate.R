suppressPackageStartupMessages(library(edgeR))

run_edgeRQLFDetRate <- function(L) {
  message("edgeRQLFDetRate")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      dge <- DGEList(L$count, group = L$condt)
      dge <- calcNormFactors(dge)
      cdr <- scale(colMeans(L$count > 0))
      design <- model.matrix(~ cdr + L$condt)
      dge <- estimateDisp(dge, design = design)
      fit <- glmQLFit(dge, design = design)
      qlf <- glmQLFTest(fit)
      tt <- topTags(qlf, n = Inf)
    })
    
    plotBCV(dge)
    plotQLDisp(fit)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
    plotSmear(qlf)
    
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  }, error = function(e) {
    "edgeRQLFDetRate results could not be calculated"
    list(session_info = session_info)
  })
}
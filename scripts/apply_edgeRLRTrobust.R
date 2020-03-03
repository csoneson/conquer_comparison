suppressPackageStartupMessages(library(edgeR))

run_edgeRLRTrobust <- function(L) {
  message("edgeRLRTrobust")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(L$count, group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~L$condt)
    dge <- estimateGLMRobustDisp(dge, design = design)
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
}
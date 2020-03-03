suppressPackageStartupMessages(library(edgeR))

run_edgeRQLF <- function(L) {
  message("edgeRQLF")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(L$count, group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~L$condt)
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
}
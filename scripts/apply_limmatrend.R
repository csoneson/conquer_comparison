suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

run_limmatrend <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(L$count, group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~L$condt)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(y, design = design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
  })
  
  hist(tt$P.Value, 50)
  hist(tt$adj.P.Val, 50)
  limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
  plotMD(fit)
  
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$P.Value,
                       padj = tt$adj.P.Val,
                       row.names = rownames(tt)))
}

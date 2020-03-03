suppressPackageStartupMessages(library(ROTS))
suppressPackageStartupMessages(library(edgeR))

run_ROTSvoom <- function(L) {
  message("ROTS, voom")
  session_info <- sessionInfo()
  timing <- system.time({
    stopifnot(all(names(L$condt) == colnames(L$count)))
    grp <- L$condt
    dge <- DGEList(counts = L$count)
    dge <- edgeR::calcNormFactors(dge)
    vm <- voom(dge, design = model.matrix(~ grp))
    rots <- ROTS(data = vm$E, groups = grp, B = 1000, K = 1000, log = TRUE, seed = 123)
  })
  
  hist(rots$pvalue, 50)
  hist(rots$FDR, 50)
  hist(rots$logfc, 50)
  print(rots$R)
  print(rots$Z)
  print(rots$k)
  print(rots$a1)
  print(rots$a2)
  
  list(session_info = session_info,
       timing = timing,
       res = rots,
       df = data.frame(pval = rots$pvalue,
                       row.names = rownames(rots$data)))
}
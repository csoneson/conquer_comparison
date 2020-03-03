suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scDD))

run_scDD <- function(L) {
  message("scDD")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      scDatList <- list()
      (groups <- unique(L$condt))
      for (i in 1:length(groups)) {
        scDatList[[paste0("G", i)]] <- as.matrix(L$count[, which(L$condt == groups[i])])
      }
      datNorm.scran <- scDD::preprocess(scDatList, 
                                        ConditionNames = names(scDatList),
                                        zero.thresh = 1, median_norm = TRUE)
      condition <- L$condt[colnames(datNorm.scran)]
      condition <- as.numeric(as.factor(condition))
      names(condition) <- colnames(datNorm.scran)
      
      SDSumExp <- SummarizedExperiment(assays = list("NormCounts" = datNorm.scran),
                                       colData = data.frame(condition))
      prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
      scd <- scDD(SDSumExp, prior_param = prior_param, testZeroes = FALSE,
                  param = BiocParallel::MulticoreParam(workers = 1), 
                  condition = "condition", min.size = 3, min.nonzero = NULL)
      res <- results(scd)
    })
    
    hist(res$nonzero.pvalue, 50)
    hist(res$nonzero.pvalue.adj, 50)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$nonzero.pvalue,
                         padj = res$nonzero.pvalue.adj,
                         row.names = rownames(res)))
  }, error = function(e) {
    "scDD results could not be calculated"
    list(session_info = session_info)
  })
}

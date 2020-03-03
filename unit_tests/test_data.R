## Test data generation and subsetting

topdir <- "/home/Shared/data/seq/conquer/comparison"

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(rjson))
source(paste0(topdir, "/scripts/prepare_mae.R"))

test_that("cells in subsets are assigned the right class", {
  gw <- getwd()
  setwd(topdir)
  for (ds in c("EMTAB2805", "EMTAB2805mock", "GSE45719", "GSE45719mock", "GSE74596", "GSE74596mock", "UsoskinGSE59739", "UsoskinGSE59739mock", "EGEUV1", "EGEUV1mock", "GSE63818-GPL16791", "GSE60749-GPL13112", "GSE60749-GPL13112mock", "GSE48968-GPL13112", "GSE48968-GPL13112mock", "GSE45719sim123", "GSE45719sim123mock", "GSE74596sim123", "GSE74596sim123mock", "GSE48968-GPL13112sim123", "GSE48968-GPL13112sim123mock")) {
    config <- fromJSON(file = paste0("config/", ds, ".json"))
    subsets <- readRDS(config$subfile)
    data <- readRDS(config$mae)
    data <- clean_mae(mae = data, groupid = config$groupid)
    pdt <- pData(data)
    for (i in config$sizes) {
      for (j in nrow(subsets$out_condition[[as.character(i)]])) {
        datasub <- subset_mae(data, subsets$keep_samples, sz = i, i = j, 
                              imposed_condition = subsets$out_condition, filt = "")
        
        ## Check that the condition in "out_condition" is the same as that given
        ## in the phenodata slot of the data set for the same sample
        expect_equal(as.character(pdt[match(subsets$keep_samples[[as.character(i)]][j, ], 
                                            rownames(pdt)), paste(config$groupid, collapse = ".")]), 
                     as.character(gsub("\\.[1-2]$", "", subsets$out_condition[[as.character(i)]][j, ])))

        ## Check that all retained samples (in "keep_samples") belong to one of
        ## the groups that are intended to be kept
        expect_equal(as.character(sort(unique(pdt[match(subsets$keep_samples[[as.character(i)]][j, ], 
                                                        rownames(pdt)), 
                                                  paste(config$groupid, collapse = ".")]))), 
                     as.character(sort(config$keepgroups)))
        
        ## Check that the condition vector of the subsetted data set matches the
        ## "out_condition"
        expect_equivalent(datasub$condt[match(subsets$keep_samples[[as.character(i)]][j, ],
                                              names(datasub$condt))],
                          subsets$out_condition[[as.character(i)]][j, ])
        
        ## Check that the cells in the subsetted objects are all in the same
        ## order
        expect_equal(colnames(datasub$count), colnames(datasub$tpm))
        expect_equal(colnames(datasub$count), names(datasub$condt))
      }
    }
  }
  setwd(gw)
})
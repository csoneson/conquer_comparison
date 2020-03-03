## Test COBRAData generation

topdir <- "/home/Shared/data/seq/conquer/comparison"

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(iCOBRA))
source(paste0(topdir, "/scripts/prepare_mae.R"))

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

test_that("COBRAData object is correctly assembled", {
  gw <- getwd()
  setwd(topdir)
  for (ds in c("EMTAB2805", "EMTAB2805mock", "GSE45719", "GSE45719mock", "GSE74596", "GSE74596mock", "UsoskinGSE59739", "UsoskinGSE59739mock", "EGEUV1", "EGEUV1mock", "GSE63818-GPL16791", "GSE60749-GPL13112", "GSE60749-GPL13112mock", "GSE48968-GPL13112", "GSE48968-GPL13112mock", "GSE45719sim123", "GSE45719sim123mock", "GSE74596sim123", "GSE74596sim123mock", "GSE48968-GPL13112sim123", "GSE48968-GPL13112sim123mock")) {
    config <- fromJSON(file = paste0("config/", ds, ".json"))
    subsets <- readRDS(config$subfile)
    data <- readRDS(config$mae)
    data <- clean_mae(mae = data, groupid = config$groupid)
    
    for (f in c("", "_TPM_1_25p")) {
      cobra <- readRDS(paste0(topdir, "/output/cobra_data/", ds, f, "_cobra.rds"))
      ngenes <- readRDS(paste0(topdir, "/output/cobra_data/", ds, f, "_nbr_called.rds"))
      
      all_methods <- unique(get_method(colnames(padj(cobra))))
      for (mth in all_methods) {
        ## Test that adjusted p-values in COBRAData object are the same as in
        ## the individual result files
        res <- readRDS(paste0(topdir, "/results/", ds, "_", mth, ".rds"))
        for (nm in names(res)) {
          if ("padj" %in% colnames(res[[nm]]$df)) {
            expect_equal(res[[nm]]$df$padj[match(rownames(padj(cobra)), rownames(res[[nm]]$df))], 
                         padj(cobra)[, nm])
          }
          datasub <- subset_mae(data, subsets$keep_samples, sz = get_nsamples(nm), 
                                i = as.numeric(as.character(get_repl(nm))), 
                                imposed_condition = subsets$out_condition, 
                                filt = gsub("^_", "", f))
          ## Test that the calculated number of tested, called and significant
          ## genes are correct
          if (nm %in% colnames(padj(cobra))) {
            ntested <- nrow(datasub$count)
            ncalled <- sum(!is.na(padj(cobra)[, nm]))
            nsign <- sum(padj(cobra)[, nm] <= 0.05, na.rm = TRUE)
            tmp <- subset(ngenes, method == mth & ncells == get_nsamples(nm) & 
                            repl == get_repl(nm) & dataset == ds & 
                            filt == gsub("^_", "", f))
            expect_equal(ntested, tmp$nbr_tested)
            expect_equal(ncalled, tmp$nbr_called)
            expect_equal(nsign, tmp$nbr_sign_adjp0.05)
          }
        }
        
        ## Test that truth is equivalent to adjusted p-values from largest
        ## sample set
        maxn <- max(as.numeric(as.character(get_nsamples(names(res)))))
        if (paste0(mth, ".", maxn, ".1") %in% colnames(padj(cobra))) {
          tmp <- padj(cobra)[, paste0(mth, ".", maxn, ".1")]
          tmp[is.na(tmp)] <- 1
          expect_equal(as.numeric(tmp <= 0.05),
                       truth(cobra)[match(rownames(padj(cobra)), rownames(truth(cobra))), 
                                    paste0(mth, ".truth")])
        }
      }
    }        
  }
  setwd(gw)
})
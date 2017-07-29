suppressPackageStartupMessages(library(scImpute))

impute_dropouts <- function(count, tpm, condt, avetxlength) {
  rnb <- paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), "_", round(runif(1) * 1e8))
  
  count <- as.matrix(count)
  totalCounts_by_cell <- colSums(count)
  totalCounts_by_cell[totalCounts_by_cell == 0] <- 1
  count <- sweep(count, MARGIN = 2, totalCounts_by_cell/10^6, FUN = "/")
  if (min(count) < 0) {
    message("smallest count cannot be negative!")
    stop()
  }
  count <- log10(count + 1.01)
  
  genenames <- rownames(count)
  cellnames <- colnames(count)
  message("estimating mixture models ...")
  scImpute:::get_mix_parameters(count = count, point = log10(1.01), 
                                path = paste0("tmp/", rnb, "parslist.rds"), 
                                ncores = 1)
  parslist <- readRDS(paste0("tmp/", rnb, "parslist.rds"))
  count_imp <- scImpute:::imputation_model1_bytype(count = count, 
                                                   labels = condt, 
                                                   point = log10(1.01), 
                                                   parslist = parslist, 
                                                   drop_thre = 0.5, method = 2, 
                                                   ncores = 1)
  count_imp <- 10 ^ count_imp - 1.01
  rownames(count_imp) <- genenames
  colnames(count_imp) <- cellnames
  
  count_imp <- sweep(count_imp, MARGIN = 2, totalCounts_by_cell/10^6, 
                     FUN = "*")

  ## Estimate TPMs
  tpm_imp <- count_imp/rowMeans(avetxlength)
  tpm_imp <- t(t(tpm_imp) / colSums(tpm_imp)) * 1e6
  
  list(count = count_imp, tpm = tpm_imp, condt = condt)
}
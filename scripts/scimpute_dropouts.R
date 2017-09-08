suppressPackageStartupMessages(library(scImpute))

## Impute dropouts with the scImpute package
scimpute_dropouts <- function(count, tpm, condt, avetxlength) {
  stopifnot(length(unique(condt)) == 2)
  rnb <- paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), "_", round(runif(1) * 1e8))
  
  count_orig <- as.matrix(count)  ## save for later comparison
  
  count <- as.matrix(count)
  totalCounts_by_cell <- colSums(count)  ## library sizes
  totalCounts_by_cell[totalCounts_by_cell == 0] <- 1
  count <- sweep(count, MARGIN = 2, totalCounts_by_cell/10^6, FUN = "/")  ## cpm
  if (min(count) < 0) {
    message("smallest count cannot be negative!")
    stop()
  }
  count <- log10(count + 1.01)  ## log10(cpm)
  genenames <- rownames(count)
  
  ## Split by condition and impute separately, as suggested in the scImpute
  ## vignette
  countlist <- list(gr1 = count[, which(condt == levels(factor(condt))[1])],
                    gr2 = count[, which(condt == levels(factor(condt))[2])])
  cellnames <- lapply(countlist, colnames)
  countimplist <- list()
  for (nm in names(countlist)) {
    message("estimating mixture models ...")
    scImpute:::get_mix_parameters(count = countlist[[nm]], point = log10(1.01), 
                                  path = paste0("tmp/", rnb, nm, "parslist.rds"), 
                                  ncores = 1)
    parslist <- readRDS(paste0("tmp/", rnb, nm, "parslist.rds"))
    countimplist[[nm]] <- scImpute:::imputation_model1(count = countlist[[nm]], 
                                                       point = log10(1.01), 
                                                       parslist = parslist, 
                                                       drop_thre = 0.5, method = 2, 
                                                       ncores = 1)
    colnames(countimplist[[nm]]) <- cellnames[[nm]]
  }
  ## Merge the two conditions
  count_imp <- do.call(cbind, countimplist)
  count_imp <- count_imp[, match(colnames(count), colnames(count_imp))]
  count_imp <- 10 ^ count_imp - 1.01
  rownames(count_imp) <- genenames

  count_imp <- sweep(count_imp, MARGIN = 2, totalCounts_by_cell/10^6, 
                     FUN = "*")

  ## Estimate TPMs
  stopifnot(all(rownames(count_imp) == rownames(avetxlength)))
  tpm_imp <- count_imp/rowMeans(avetxlength)
  tpm_imp <- t(t(tpm_imp) / colSums(tpm_imp)) * 1e6
  
  ## Tabulate number of imputed values
  stopifnot(all(colnames(count_imp) == colnames(count_orig)))
  stopifnot(all(rownames(count_imp) == rownames(count_orig)))
  nimp <- data.frame(gene = rownames(count_imp), 
                     nbr_increased = rowSums(count_imp > (count_orig + 1e-6)),
                     nbr_decreased = rowSums(count_imp < (count_orig - 1e-6))) %>%
    dplyr::mutate(nbr_unchanged = ncol(count_imp) - nbr_increased - nbr_decreased)
  rownames(nimp) <- nimp$gene
  
  list(count = count_imp, tpm = tpm_imp, condt = condt, nimp = nimp)
}
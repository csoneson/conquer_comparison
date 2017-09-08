suppressPackageStartupMessages(library(DrImpute))
suppressPackageStartupMessages(library(dplyr))

## Impute dropouts with the DrImpute package
drimpute_dropouts <- function(count, tpm, condt, avetxlength) {
  ## Define ks to use
  n <- ncol(count)/2
  ncl <- 6
  ks <- unique(round(2:n))
  if (length(ks) <= ncl) {
    ks <- ks
  } else if (length(ks) <= (3 * ncl)) {
    ks <- ks[ks <= (4 * ncl)]
    ks <- ks[((length(ks) - ncl)/2 + 1):(length(ks) - (n - ncl)/2)]
  } else {
    ks <- ks[ks <= (4 * ncl)]
    ks <- ks[seq(1, length(ks), 2)]
    ks <- ks[((length(ks) - ncl)/2 + 1):(length(ks) - (length(ks) - ncl)/2)]
  }
  ks
  
  count_orig <- as.matrix(count)  ## save for later comparison
  
  ## Normalize counts
  count <- as.matrix(count)
  sf <- colSums(count)
  sf[sf == 0] <- 1
  sf <- sf/exp(mean(log(sf)))
  lognormcount <- log(sweep(count, 2, sf, "/") + 1)

  lognormcount_imp <- 
    DrImpute(lognormcount, ks = ks, dists = c("spearman", "pearson"), 
             method = "mean", cls = NULL, seed = 123, zerop = 0)
  
  ## Go back to original count scale
  normcount_imp <- exp(lognormcount_imp) - 1
  count_imp <- sweep(normcount_imp, 2, sf, "*")
  colnames(count_imp) <- colnames(count_orig)
  
  ## Estimate TPMs
  stopifnot(!is.null(colnames(count_imp)))
  stopifnot(!is.null(rownames(count_imp)))
  stopifnot(all(rownames(count_imp) == rownames(avetxlength)))
  tpm_imp <- count_imp/rowMeans(avetxlength)
  tpm_imp <- t(t(tpm_imp) / colSums(tpm_imp)) * 1e6
  
  ## Tabulate number of imputed values
  stopifnot(!is.null(colnames(tpm_imp)))
  stopifnot(!is.null(rownames(tpm_imp)))
  stopifnot(all(colnames(count_imp) == colnames(count_orig)))
  stopifnot(all(rownames(count_imp) == rownames(count_orig)))
  stopifnot(all(colnames(tpm_imp) == colnames(count_orig)))
  stopifnot(all(rownames(tpm_imp) == rownames(count_orig)))
  nimp <- data.frame(gene = rownames(count_imp), 
                     nbr_increased = rowSums(count_imp > (count_orig + 1e-6)),
                     nbr_decreased = rowSums(count_imp < (count_orig - 1e-6))) %>%
    dplyr::mutate(nbr_unchanged = ncol(count_imp) - nbr_increased - nbr_decreased)
  rownames(nimp) <- nimp$gene
  
  list(count = count_imp, tpm = tpm_imp, condt = condt, nimp = nimp)
}
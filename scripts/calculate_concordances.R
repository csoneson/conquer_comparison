args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

source("scripts/concordance_functions.R")

print(dataset)  ## Data set
print(filt)  ## Filtering
print(cobradir)  ## Directory where to look for cobra object (output from prepare_cobra_for_evaluation.R)
print(outdir) ## Directory where to save the calculated concordances

if (filt == "") {
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

## Set maximal number of top-ranked genes to consider
maxrank <- 1000

## ----------------------------- Help functions ----------------------------- ##
get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

## ---------------------------- Prepare data -------------------------------- ##

cobra <- readRDS(paste0(cobradir, "/", dataset, exts, "_cobra.rds"))
## Set p-values and adjusted p-values for untested genes and genes with NA
## values to 1, so that they will be ranked last
pval(cobra)[is.na(pval(cobra))] <- 1
padj(cobra)[is.na(padj(cobra))] <- 1
## Similarly, set score for untested genes and genes with NA values to a value 
## below the smallest observed score and change the sign so that the smallest
## (most negative) score corresponds to the most significant gene
if (nrow(iCOBRA::score(cobra)) > 0) {
  iCOBRA::score(cobra)[is.na(iCOBRA::score(cobra))] <- 
    min(iCOBRA::score(cobra), na.rm = TRUE) - 1
  iCOBRA::score(cobra) <- -iCOBRA::score(cobra)
}

summary_data <- list()

## Define the set of values that will be used to rank the genes
pconc <- pval(cobra)
## For methods not returning p-values, use scores or adjusted p-values
addm <- setdiff(colnames(iCOBRA::score(cobra)), colnames(pconc))
if (length(addm) > 0) {
  pconc <- dplyr::full_join(data.frame(gene = rownames(pconc), pconc),
                            data.frame(gene = rownames(iCOBRA::score(cobra)), 
                                       iCOBRA::score(cobra)[, addm, drop = FALSE]))
  rownames(pconc) <- pconc$gene
  pconc$gene <- NULL
}
addm <- setdiff(colnames(padj(cobra)), colnames(pconc))
if (length(addm) > 0) {
  pconc <- dplyr::full_join(data.frame(gene = rownames(pconc), pconc),
                            data.frame(gene = rownames(padj(cobra)), 
                                       padj(cobra)[, addm, drop = FALSE]))
  rownames(pconc) <- pconc$gene
  pconc$gene <- NULL
}
## Find ordering of each column and convert to matrix
for (i in colnames(pconc)) {
  pconc[, i] <- order(pconc[, i])
}
pconc <- as.matrix(pconc)

all_methods <- unique(get_method(colnames(pconc)))
all_nsamples <- unique(get_nsamples(colnames(pconc)))
all_repl <- unique(get_repl(colnames(pconc)))

## ----------------------- Calculate concordances --------------------------- ##
## Calculate the number of shared occurrences among top-k genes across a number
## of different scenarios. Also calculate concordance values (number of genes
## shared by *all* columns for given k) and calculate partial AUC for the
## concordance curve (k vs concordance). The AUCs is a scaled AUC, obtained by
## dividing with the maximum area = k^2/2

## Across all instances (all sample sizes, all replicates), for given method
nbrshared <- do.call(rbind, lapply(all_methods, function(mth) {
  tmp <- calculate_nbr_occurrences(mtx = pconc[, which(get_method(colnames(pconc)) == mth), 
                                               drop = FALSE], 
                                   maxrank = maxrank)
  if (!is.null(tmp)) {
    tmp$method <- mth
    tmp
  } else {
    NULL
  }
})) %>% dplyr::mutate(dataset = dataset, filt = filt)
concvals <- nbrshared %>% 
  dplyr::filter(nbr_occ == nbr_cols) %>%
  dplyr::group_by(method) %>%
  calc_auc()
summary_data$concordance_fullds_bymethod <- concvals
summary_data$nbrshared_fullds_bymethod <- nbrshared

## Between pairwise instances with the same sample size, for given method
nbrshared_pairwise <- do.call(rbind, lapply(all_methods, function(mth) {
  do.call(rbind, lapply(all_nsamples, function(i) {
    tmp <- pconc[, intersect(which(get_method(colnames(pconc)) == mth),
                             which(get_nsamples(colnames(pconc)) == i)), drop = FALSE]
    if (ncol(tmp) > 1) {
      nsr <- NULL
      for (j1 in 1:(ncol(tmp) - 1)) {
        for (j2 in (j1 + 1):(ncol(tmp))) {
          cv <- calculate_nbr_occurrences(mtx = tmp[, c(j1, j2)], maxrank = maxrank)
          cv$ncells <- i
          cv$replicate1 <- get_repl(colnames(tmp)[j1])
          cv$replicate2 <- get_repl(colnames(tmp)[j2])
          nsr <- rbind(nsr, cv)
        }
      }
      nsr$method <- mth
      nsr
    } else {
      NULL
    }
  }))
})) %>% dplyr::mutate(dataset = dataset, filt = filt)
concvals_pairwise <- nbrshared_pairwise %>% 
  dplyr::filter(nbr_occ == nbr_cols) %>%
  dplyr::group_by(method, ncells, replicate1, replicate2) %>%
  calc_auc()
summary_data$concordance_pairwise_bymethod <- concvals_pairwise
summary_data$nbrshared_pairwise_bymethod <- nbrshared_pairwise

## Between pairs of methods, for a given data set instance (fixed sample size, replicate)
nbrshared_btwmth <- do.call(rbind, lapply(all_nsamples, function(ss) {
  do.call(rbind, lapply(all_repl, function(i) {
    tmp <- pconc[, intersect(which(get_repl(colnames(pconc)) == i),
                             which(get_nsamples(colnames(pconc)) == ss)), drop = FALSE]
    if (ncol(tmp) > 1) {
      nsr <- NULL
      for (j1 in 1:(ncol(tmp) - 1)) {
        for (j2 in (j1 + 1):(ncol(tmp))) {
          cv <- calculate_nbr_occurrences(mtx = tmp[, c(j1, j2)], maxrank = maxrank)
          cv$method1 <- get_method(colnames(tmp)[j1])
          cv$method2 <- get_method(colnames(tmp)[j2])
          nsr <- rbind(nsr, cv)
        }
      }
      nsr$ncells <- ss
      nsr$repl <- i
      nsr
    } else {
      NULL
    }
  }))
})) %>% dplyr::mutate(dataset = dataset, filt = filt)
concvals_btwmth <- nbrshared_btwmth %>% 
  dplyr::filter(nbr_occ == nbr_cols) %>%
  dplyr::group_by(method1, method2, ncells, repl) %>%
  calc_auc()
summary_data$concordance_betweenmethods_pairwise <- concvals_btwmth
summary_data$nbrshared_betweenmethods_pairwise <- nbrshared_btwmth

## -------------------------- Save output ----------------------------------- ##
saveRDS(summary_data, file = paste0(outdir, "/", dataset, exts, "_concordances.rds"))

sessionInfo()
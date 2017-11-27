## Test concordance calculations with small example

suppressPackageStartupMessages(library(iCOBRA))

topdir <- "/home/Shared/data/seq/conquer/comparison"

source(paste0(topdir, "/scripts/concordance_functions.R"))

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

test_that("concordance calculations are correct on small example", {
  R <- matrix(c(1, 5, 2, 2, 8, 1, 3, 3, 7, 4, 9, 8, 5, 6, 3), byrow = TRUE, nrow = 5)
  nbr <- calculate_nbr_occurrences(mtx = R, maxrank = 5)
  expect_equal(nbr$nbr_occ, rep(1:3, 5))
  expect_equal(nbr$k, rep(1:5, each = 3))
  expect_equal(nbr$nbr_genes, c(3, 0, 0, 2, 2, 0, 3, 3, 0, 4, 4, 0, 4, 4, 1))
  expect_equal(nbr$nbr_cols, rep(3, 15))
  
  auc <- calc_auc(nbr %>% dplyr::filter(nbr_occ == nbr_cols))
  expect_equal(auc$nbr_occ, rep(3, 5))
  expect_equal(auc$k, 1:5)
  expect_equal(auc$nbr_genes, c(0, 0, 0, 0, 1))
  expect_equal(auc$AUC, c(0, 0, 0, 0, 0.5))
  expect_equal(auc$AUCs, c(0, 0, 0, 0, 0.04))
  
  R <- cbind(1:15, 1:15)
  nbr <- calculate_nbr_occurrences(mtx = R, maxrank = 15)
  auc <- calc_auc(nbr %>% dplyr::filter(nbr_occ == nbr_cols))
  expect_equal(auc$AUCs, rep(1, 15))
})

test_that("concordance calculations are correct on real data", {
  ds <- "GSE74596"
  exts <- "_TPM_1_25p"
  cobra <- readRDS(paste0(topdir, "/output/cobra_data/", ds, exts, "_cobra.rds"))
  concordances <- readRDS(paste0(topdir, "/output/concordances/", ds, exts, "_concordances.rds"))

  pval(cobra)[is.na(pval(cobra))] <- 1
  padj(cobra)[is.na(padj(cobra))] <- 1
  if (nrow(iCOBRA::score(cobra)) > 0) {
    iCOBRA::score(cobra)[is.na(iCOBRA::score(cobra))] <- 
      min(iCOBRA::score(cobra), na.rm = TRUE) - 1
    iCOBRA::score(cobra) <- -iCOBRA::score(cobra)
  }
  
  ## Define the set of values that will be used to rank the genes
  pv <- pval(cobra)
  ## For methods not returning p-values, use scores or adjusted p-values
  addm <- setdiff(colnames(iCOBRA::score(cobra)), colnames(pv))
  if (length(addm) > 0) {
    pv <- dplyr::full_join(data.frame(gene = rownames(pv), pv),
                              data.frame(gene = rownames(iCOBRA::score(cobra)), 
                                         iCOBRA::score(cobra)[, addm]))
    rownames(pv) <- pv$gene
    pv$gene <- NULL
  }
  addm <- setdiff(colnames(padj(cobra)), colnames(pv))
  if (length(addm) > 0) {
    pv <- dplyr::full_join(data.frame(gene = rownames(pv), pv),
                              data.frame(gene = rownames(padj(cobra)), 
                                         padj(cobra)[, addm]))
    rownames(pv) <- pv$gene
    pv$gene <- NULL
  }
  ## Find ordering of each column and convert to matrix
  for (i in colnames(pv)) {
    pv[, i] <- order(pv[, i])
  }
  pv <- as.matrix(pv)
  
  maxn <- 100
  
  ## Single method, all instances
  mth <- "DESeq2"
  tmp <- pv[1:maxn, which(get_method(colnames(pv)) == paste0(mth, exts))]
  expect_equal(ncol(tmp), 16)
  
  conc <- data.frame(t(sapply(1:maxn, function(i) {
    p <- sort(tmp[1:i, ])
    p <- sum(table(p) == ncol(tmp))
    c(k = i, p = p)
  })), stringsAsFactors = FALSE)
  
  trueconc <- concordances$concordance_fullds_bymethod %>% 
    dplyr::filter(method == paste0(mth, exts), dataset == ds)
  expect_equal(trueconc$k[1:maxn], conc$k)
  expect_equal(trueconc$nbr_genes[1:maxn], conc$p)
  
  conc_auc <- caTools::trapz(c(0, conc$k, conc$k[length(conc$k)]), 
                             c(0, conc$p, 0))/(maxn^2/2)
  expect_equal(trueconc$AUCs[trueconc$k == maxn], conc_auc)
  
  ## Single method, given pair of instances
  mth <- "DESeq2"
  n_cells <- 12
  inst1 <- 1
  inst2 <- 2
  tmp <- pv[1:maxn, which(get_method(colnames(pv)) == paste0(mth, exts) & 
                            get_nsamples(colnames(pv)) == n_cells & 
                            get_repl(colnames(pv)) %in% c(inst1, inst2))]
  expect_equal(ncol(tmp), 2)
  
  conc <- data.frame(t(sapply(1:maxn, function(i) {
    p <- sort(tmp[1:i, ])
    p <- sum(table(p) == ncol(tmp))
    c(k = i, p = p)
  })), stringsAsFactors = FALSE)
  
  trueconc <- concordances$concordance_pairwise_bymethod %>% 
    dplyr::filter(method == paste0(mth, exts)) %>%
    dplyr::filter(dataset == ds) %>%
    dplyr::filter(ncells == n_cells) %>%
    dplyr::filter(replicate1 == inst1) %>%
    dplyr::filter(replicate2 == inst2)
  expect_equal(trueconc$k[1:maxn], conc$k)
  expect_equal(trueconc$nbr_genes[1:maxn], conc$p)
  
  conc_auc <- caTools::trapz(c(0, conc$k, conc$k[length(conc$k)]), 
                             c(0, conc$p, 0))/(maxn^2/2)
  expect_equal(trueconc$AUCs[trueconc$k == maxn], conc_auc)
  
  ## Pair of methods, given instance
  mth1 <- "DESeq2"
  mth2 <- "Wilcoxon"
  n_cells <- 12
  inst1 <- 1
  tmp <- pv[1:maxn, which(get_method(colnames(pv)) %in% paste0(c(mth1, mth2), exts) & 
                            get_nsamples(colnames(pv)) == n_cells & 
                            get_repl(colnames(pv)) == inst1)]
  expect_equal(ncol(tmp), 2)
  
  conc <- data.frame(t(sapply(1:maxn, function(i) {
    p <- sort(tmp[1:i, ])
    p <- sum(table(p) == ncol(tmp))
    c(k = i, p = p)
  })), stringsAsFactors = FALSE)
  
  trueconc <- concordances$concordance_betweenmethods_pairwise %>% 
    dplyr::filter(dataset == ds) %>%
    dplyr::filter(ncells == n_cells) %>%
    dplyr::filter(repl == inst1) %>%
    dplyr::filter(method1 == paste0(mth2, exts)) %>%
    dplyr::filter(method2 == paste0(mth1, exts))
  expect_equal(trueconc$k[1:maxn], conc$k)
  expect_equal(trueconc$nbr_genes[1:maxn], conc$p)
  
  conc_auc <- caTools::trapz(c(0, conc$k, conc$k[length(conc$k)]), 
                             c(0, conc$p, 0))/(maxn^2/2)
  expect_equal(trueconc$AUCs[trueconc$k == maxn], conc_auc)
  
})
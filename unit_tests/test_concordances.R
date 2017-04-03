## Test concordance calculations with small example

topdir <- "/home/Shared/data/seq/conquer/comparison"

source(paste0(topdir, "/scripts/concordance_functions.R"))

test_that("concordance calculations are correct", {
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
})
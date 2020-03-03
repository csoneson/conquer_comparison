## For each k in 1:maxrank, count the number of genes occurring each number of 
## times. The output is a data frame with three columns: nbr_occ, nbr_genes, k. 
## For a given row, interpret as follows: among the top-k genes from each 
## method, nbr_genes occur exactly nbr_occ times. 
calculate_nbr_occurrences <- function(mtx, maxrank) {
  maxrank <- min(maxrank, nrow(mtx))
  if (ncol(mtx) > 1) {
    M <- matrix(0, max(mtx[1:maxrank, ]), maxrank)
    for (i in 1:ncol(mtx)) {
      M[cbind(mtx[1:maxrank, i], 1:maxrank)] <- M[cbind(mtx[1:maxrank, i], 1:maxrank)] + 1
    }
    M <- M[rowSums(M) != 0, ]
    M <- t(apply(M, 1, cumsum))
    M2 <- matrix(0, ncol(mtx), maxrank)
    for (i in 1:nrow(M)) {
      M2[cbind(M[i, ], 1:ncol(M))] <- M2[cbind(M[i, ], 1:ncol(M))] + 1
    }
    M2 <- M2 %>% reshape2::melt() %>%
      dplyr::rename(k = Var2, nbr_genes = value, nbr_occ = Var1)
    M2 %>% dplyr::arrange(k, nbr_occ) %>%
      dplyr::mutate(nbr_cols = ncol(mtx))
  } else {
    NULL
  }
}

## Calculate partial (cumulative) AUCs. 
## Assumes that x variable = k, y variable = nbr_genes
calc_auc <- function(x) {
  x %>% dplyr::mutate(dx = c(k[1], diff(k)),
                      dy = c(nbr_genes[1], diff(nbr_genes)),
                      ys = c(0, nbr_genes[-length(nbr_genes)])) %>%
    dplyr::mutate(AUC = cumsum(dx * dy/2 + dx * ys)) %>%
    dplyr::mutate(AUCs = AUC/(k^2/2))
}

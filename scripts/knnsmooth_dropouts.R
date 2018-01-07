## Impute dropouts using the knn-smooth algorithm (code from https://github.com/yanailab/knn-smoothing)

## The following code is from the above GitHub repository, and was copied on December 14, 2017
# Author: Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017 New York University

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(magrittr))

freeman_tukey_transform <- function(mat) {
  sqrt(mat) + sqrt(mat + 1)
}

calculate_distances <- function(mat) {
  # mat: gene by sample
  # normalize to median transcript count
  num_transcripts <- Matrix::colSums(mat)
  size_factor <- median(num_transcripts, na.rm = T) / num_transcripts
  
  mat_norm <- t(t(mat) * size_factor)
  # apply freeman-tukey transform
  mat_FTT <- freeman_tukey_transform(mat_norm)
  # calculate all pairwise distances using the Euclidean metric
  mat_D <- dist(t(mat_FTT), method = "euclidean",
                upper = T, diag = T)
  return(as.matrix(mat_D))
}


#' KNN-smoothing on UMI-filtered single-cell RNA-seq data
#'
#' @param mat A numeric matrix with gene names on rows and cell names on columns.
#' @param k Number of nearest neighbours to aggregate.
#'
#' @return A smoothed numeric matrix.
#' @export
#'
#' @examples
#' X <- matrix(abs(sin(seq(from=1, to=1000, length.out = 1000))),
#' nrow = 25, byrow = T)
#' y <- rep(1:4, each=10)
#' dim(X)
#' colnames(X) <- as.character(paste0('s', seq_len(ncol(X))))
#' rownames(X) <- as.character(paste0('g', seq_len(nrow(X))))
#' S <- knn_smoothing(X, k=5)
#' plot(X[1, ], X[3, ], col=factor(y), main='original')
#' plot(S[1, ], S[3, ], col=factor(y), main='smoothed')
knn_smoothing <- function(mat, k=5) {
  cname <- colnames(mat)
  gname <- rownames(mat)
  num_powers <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_powers)){
    k_step <- min(2^p - 1, k)
    message(paste0('Step ', p, '/', num_powers, ':',
                   'Smoothing using k=', k_step))
    D <- calculate_distances(S)
    S <- sapply(cname, function(cn){
      closest_id <- D[cn, ] %>% sort(.) %>%
        head(., k_step+1) %>%
        names(.)
      closest_mat <- mat[gname, closest_id] %>%
        matrix(., nrow=length(gname), byrow = F)
      rownames(closest_mat) <- gname
      colnames(closest_mat) <- closest_id
      return(Matrix::rowSums(closest_mat))
    })
  }
  S
}

knnsmooth_dropouts <- function(count, tpm, condt, avetxlength) {
  ## Determine k
  k <- min(length(condt)/4, 15)
  
  ## Smooth counts
  count_imp <- knn_smoothing(count, k = k)
  
  ## Estimate TPMs
  stopifnot(!is.null(colnames(count_imp)))
  stopifnot(!is.null(rownames(count_imp)))
  stopifnot(all(rownames(count_imp) == rownames(avetxlength)))
  tpm_imp <- count_imp/rowMeans(avetxlength)
  tpm_imp <- t(t(tpm_imp) / colSums(tpm_imp)) * 1e6
  
  ## Tabulate number of imputed values
  stopifnot(!is.null(colnames(tpm_imp)))
  stopifnot(!is.null(rownames(tpm_imp)))
  stopifnot(all(colnames(count_imp) == colnames(count)))
  stopifnot(all(rownames(count_imp) == rownames(count)))
  stopifnot(all(colnames(tpm_imp) == colnames(count)))
  stopifnot(all(rownames(tpm_imp) == rownames(count)))
  nimp <- data.frame(gene = rownames(count_imp), 
                     nbr_increased = rowSums(count_imp > (count + 1e-6)),
                     nbr_decreased = rowSums(count_imp < (count - 1e-6))) %>%
    dplyr::mutate(nbr_unchanged = ncol(count_imp) - nbr_increased - nbr_decreased)
  rownames(nimp) <- nimp$gene
  
  list(count = count_imp, tpm = tpm_imp, condt = condt, nimp = nimp)
}
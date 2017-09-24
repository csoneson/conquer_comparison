## Prepare a MultiAssayExperiment for the 10x Genomics data

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))

## This function is a simplified variant of the read10xResults() function in the
## scater package 
## (https://github.com/davismcc/scater/blob/master/R/10ximport-wrapper.R)
read10X <- function(data_dir) {
  data_mat <- Matrix::readMM(file.path(data_dir, "matrix.mtx"))
  keep_barcode <- Matrix::colSums(data_mat) >= 1000
  data_mat <- data_mat[, keep_barcode]
  cell.names <- data.table::fread(file.path(data_dir, "barcodes.tsv"), header = FALSE)[[1]]
  cell.names <- gsub("-.*$", "", cell.names)
  cell.names <- cell.names[keep_barcode]
  gene_info <- data.table::fread(file.path(data_dir, "genes.tsv"), header = FALSE)
  colnames(data_mat) <- cell.names
  rownames(data_mat) <- gene_info[[1]]
  as.matrix(data_mat)
}

## Read count matrices
monocytes <- read10X(data_dir = "data/10xGenomics/cd14_monocytes_matrices_mex/hg19")
cytotoxict <- read10X(data_dir = "data/10xGenomics/cytotoxic_t_matrices_mex/hg19")

## Subsample 400 cells of each type
set.seed(123)
monocytes <- monocytes[, sample(seq_len(ncol(monocytes)), 400, replace = FALSE)]
cytotoxict <- cytotoxict[, sample(seq_len(ncol(cytotoxict)), 400, replace = FALSE)]

## Modify colnames to make them unique
colnames(monocytes) <- paste0("monocytes", colnames(monocytes))
colnames(cytotoxict) <- paste0("cytotoxt", colnames(cytotoxict))

## Merge count matrices and generate metadata
counts <- cbind(monocytes, cytotoxict)
meta <- data.frame(id = colnames(counts), 
                   group = rep(c("monocyte", "cytotoxict"), 
                               c(ncol(monocytes), ncol(cytotoxict))),
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$id

## Generate MultiAssayExperiment
generse <- SummarizedExperiment(
  assays = list(count = counts,
                count_lstpm = counts,
                TPM = sweep(counts, 2, colSums(counts)/1e6, "/"))
  )

stopifnot(all(colnames(generse) == rownames(meta)))
mae <- MultiAssayExperiment(experiments = list(gene = generse),
                            pData = droplevels(meta))
mae@metadata <- list(organism = "Homo sapiens")

summary(colSums(assays(experiments(mae)[["gene"]])[["count_lstpm"]]))
summary(colSums(assays(experiments(mae)[["gene"]])[["TPM"]]))
table(pData(mae)$group)

saveRDS(mae, file = "data/10XMonoCytoT.rds")
saveRDS(mae, file = "data/10XMonoCytoTmock.rds")

sessionInfo()
date()

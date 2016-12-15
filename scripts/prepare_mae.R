clean_mae <- function(mae, groupid) {
  ## Exclude ERCC spike-ins and rescale TPMs if needed
  mae2 <- subsetByRow(mae, grep("^ERCC-", unique(unlist(rownames(mae))), 
                                invert = TRUE, value = TRUE))
  for (m in names(experiments(mae2))) {
    assays(experiments(mae2)[[m]])[["TPM"]] <- 
      sweep(assays(experiments(mae2)[[m]])[["TPM"]], 
            2, colSums(assays(experiments(mae2)[[m]])[["TPM"]]), "/") * 1e6
  }
  
  if (length(groupid) > 1) {
    Biobase::pData(mae2)[, paste(groupid, collapse = ".")] <- 
      as.character(interaction(as.data.frame(Biobase::pData(mae2)[, groupid])))
    groupid <- paste(groupid, collapse = ".")
  }
  
  mae2
}
  
subset_mae <- function(mae, keep_samples, sz, i, imposed_condition) {
  s <- keep_samples[[as.character(sz)]][i, ]
  
  ## Subset and filter data matrices
  count <- assays(experiments(mae)[["gene"]])[["count_lstpm"]][, s]
  tpm <- assays(experiments(mae)[["gene"]])[["TPM"]][, s]
  if (!is.null(imposed_condition)) {
    condt <- structure(imposed_condition[[as.character(sz)]][i, ],
                       names = rownames(Biobase::pData(mae)[s, ]))
  } else {
    condt <- structure(as.character(Biobase::pData(mae)[s, groupid]),
                       names = rownames(Biobase::pData(mae)[s, ]))
  }
  count <- count[rowSums(count) > 0, ]
  tpm <- tpm[rowSums(tpm) > 0, ]
  stopifnot(all(names(condt) == colnames(count)))
  stopifnot(all(names(condt) == colnames(tpm)))
  stopifnot(length(unique(condt)) == 2)
  
  summary(colSums(count))
  summary(rowSums(count))
  summary(rowSums(tpm))
  
  list(count = count, tpm = tpm, condt = condt)
}

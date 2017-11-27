source("scripts/imputation_scripts.R")

clean_mae <- function(mae, groupid) {
  mae@sampleMap$assay <- factor(mae@sampleMap$assay)
  mae <- updateObject(mae)
  ## mae@sampleMap$assay <- factor(mae@sampleMap$assay)
  ## Exclude ERCC spike-ins and rescale TPMs if needed
  mae2 <- subsetByRow(mae, grep("^ERCC-", unique(unlist(rownames(mae))), 
                                invert = TRUE, value = TRUE))
  for (m in names(experiments(mae2))) {
    assays(experiments(mae2)[[m]])[["TPM"]] <- 
      sweep(assays(experiments(mae2)[[m]])[["TPM"]], 
            2, colSums(assays(experiments(mae2)[[m]])[["TPM"]]), "/") * 1e6
  }
  
  if (length(groupid) > 1) {
    if (paste0(R.Version()$major, ".", R.Version()$minor) < "3.4") {
      Biobase::pData(mae2)[, paste(groupid, collapse = ".")] <- 
        as.character(interaction(as.data.frame(Biobase::pData(mae2)[, groupid])))
      groupid <- paste(groupid, collapse = ".")
    } else {
      colData(mae2)[, paste(groupid, collapse = ".")] <- 
        as.character(interaction(as.data.frame(colData(mae2)[, groupid])))
      groupid <- paste(groupid, collapse = ".")
    }
  }
  
  mae2
}
  
subset_mae <- function(mae, keep_samples, sz, i, imposed_condition, filt, 
                       groupid = NULL, impute = NULL) {
  mae@sampleMap$assay <- factor(mae@sampleMap$assay)
  mae <- updateObject(mae)
  ## mae@sampleMap$assay <- factor(mae@sampleMap$assay)
  s <- keep_samples[[as.character(sz)]][i, ]
  
  ## Subset and filter data matrices
  count <- assays(experiments(mae)[["gene"]])[["count_lstpm"]][, s]
  tpm <- assays(experiments(mae)[["gene"]])[["TPM"]][, s]
  ## Scaling factor from TPM to count_lstpm was calculated from avetxlength of all samples
  avetxlength = assays(experiments(mae)[["gene"]])[["avetxlength"]]
  if (!is.null(imposed_condition)) {
    if (paste0(R.Version()$major, ".", R.Version()$minor) < "3.4") {
      condt <- structure(imposed_condition[[as.character(sz)]][i, ],
                         names = rownames(Biobase::pData(mae)[s, ]))
    } else {
      condt <- structure(imposed_condition[[as.character(sz)]][i, ],
                         names = rownames(colData(mae)[s, ]))
    }
  } else {
    if (is.null(groupid)) stop("Must provide groupid")
    if (paste0(R.Version()$major, ".", R.Version()$minor) < "3.4") {
      condt <- structure(as.character(Biobase::pData(mae)[s, groupid]),
                         names = rownames(Biobase::pData(mae)[s, ]))
    } else {
      condt <- structure(as.character(colData(mae)[s, groupid]),
                         names = rownames(colData(mae)[s, ]))
    }
  }
  
  if (!is.null(impute) && impute != "no" && !is.na(impute)) {
    imputed <- impute_dropouts(count = count, tpm = tpm, condt = condt, 
                               avetxlength = avetxlength, 
                               imputationmethod = impute)
    count <- imputed$count
    tpm <- imputed$tpm
    condt <- imputed$condt
    nimp <- imputed$nimp  ## number of imputed values
  } else {
    nimp <- NULL
  }
  
  if (filt == "") {
    count <- count[rowSums(count) > 0, ]
    tpm <- tpm[rowSums(tpm) > 0, ]
  } else {
    filt <- strsplit(filt, "_")[[1]]
    if (substr(filt[3], nchar(filt[3]), nchar(filt[3])) == "p") {
      (nbr <- as.numeric(gsub("p", "", filt[3]))/100 * ncol(count))
    } else {
      (nbr <- as.numeric(filt[3]))
    }
    if (filt[1] == "count") {
      keep_rows <- rownames(count)[which(rowSums(count > as.numeric(filt[2])) 
                                         > nbr)]
    } else if (filt[1] == "TPM") {
      keep_rows <- rownames(tpm)[which(rowSums(tpm > as.numeric(filt[2])) 
                                       > nbr)]
    } else {
      stop("First element of filt must be 'count' or 'TPM'.")
    }
    count <- count[match(keep_rows, rownames(count)), ]
    tpm <- tpm[match(keep_rows, rownames(tpm)), ]
  }
  stopifnot(all(names(condt) == colnames(count)))
  stopifnot(all(names(condt) == colnames(tpm)))
  stopifnot(length(unique(condt)) == 2)
  
  summary(colSums(count))
  summary(rowSums(count))
  summary(rowSums(tpm))
  
  if (!is.null(nimp)) {
    nimp <- nimp[match(rownames(count), rownames(nimp)), ]
  }
  
  list(count = count, tpm = tpm, condt = condt, nimp = nimp)
}

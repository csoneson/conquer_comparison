clean_mae <- function(mae, groupid) {
  mae <- updateObject(mae)
  mae@sampleMap$assay <- factor(mae@sampleMap$assay)
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
  
subset_mae <- function(mae, keep_samples, sz, i, imposed_condition, filt, 
                       groupid = NULL, impute = NULL) {
  mae <- updateObject(mae)
  mae@sampleMap$assay <- factor(mae@sampleMap$assay)
  s <- keep_samples[[as.character(sz)]][i, ]
  
  ## Subset and filter data matrices
  count <- assays(experiments(mae)[["gene"]])[["count_lstpm"]][, s]
  tpm <- assays(experiments(mae)[["gene"]])[["TPM"]][, s]
  avetxlength = assays(experiments(mae)[["gene"]])[["avetxlength"]][, s]
  if (!is.null(imposed_condition)) {
    condt <- structure(imposed_condition[[as.character(sz)]][i, ],
                       names = rownames(Biobase::pData(mae)[s, ]))
  } else {
    if (is.null(groupid)) stop("Must provide groupid")
    condt <- structure(as.character(Biobase::pData(mae)[s, groupid]),
                       names = rownames(Biobase::pData(mae)[s, ]))
  }
  
  if (!is.null(impute) && impute == "yes") {
    imputed <- impute_dropouts(count = count, tpm = tpm, condt = condt, 
                               avetxlength = avetxlength)
    count <- imputed$count
    tpm <- imputed$count
    condt <- imputed$condt
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
  
  list(count = count, tpm = tpm, condt = condt)
}

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
source("scripts/prepare_mae.R")
source("scripts/calculate_gene_characteristics.R")

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

demethods <- strsplit(demethods, ",")[[1]]

print(demethods) ## DE methods to include
print(resdir)  ## Directory from which to fetch the DE results
print(dataset)  ## Data set
print(config_file)  ## Configuration file
print(filt)  ## Filtering
print(distrdir)  ## Directory from which to fetch the distribution fit results
print(outdir)  ## Directory where the output will be written

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

## Create iCOBRA object from the result files for the different methods
(resfiles <- paste0(resdir, "/", dataset, "_", demethods, exts, ".rds"))
stopifnot(all(file.exists(resfiles)))
cobra <- NULL
timings <- list()
runstatus <- list()
for (rfn in resfiles) {   ## for each DE method
  rf <- readRDS(rfn)
  for (nm in names(rf)) {   ## for each data set instance (nm = method.ncells.repl)
    print(names(rf[[nm]]))
    
    ## Get stored timing information and results
    timings[[nm]] <- rf[[nm]]$timing
    df <- rf[[nm]]$df
    
    if (!is.null(df)) {
      runstatus[[nm]] <- "success"
    } else {
      runstatus[[nm]] <- "failure"
    }
    
    ## Add to iCOBRA object
    if ("pval" %in% colnames(df)) {
      cobra <- COBRAData(pval = setNames(data.frame(mt = df$pval,
                                                    row.names = rownames(df)), nm),
                         object_to_extend = cobra)
    }
    if ("padj" %in% colnames(df)) {
      cobra <- COBRAData(padj = setNames(data.frame(mt = df$padj,
                                                    row.names = rownames(df)), nm),
                         object_to_extend = cobra)
    }
  }
}

cobra <- calculate_adjp(cobra)

## ------------ Add gene characteristics to the COBRA object ---------------- ##

config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)
subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition

## Calculate gene characteristics for each data subset
sizes <- names(keep_samples)
truth <- list()
for (sz in sizes) {   ## for each sample size
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {   ## for each replicate
    message(sz, ".", i)
    L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i, 
                    imposed_condition = imposed_condition, filt = filt,
                    impute = config$impute)
    chars <- calculate_gene_characteristics(L)
    characs <- chars$characs
    colnames(characs)[colnames(characs) != "gene"] <- 
      paste0(colnames(characs)[colnames(characs) != "gene"], ".", sz, ".", i)
    truth[[paste0(sz, ".", i)]] <- characs
  }
}
truth <- Reduce(function(...) dplyr::full_join(..., by = "gene"), truth)

## ------------ Add information about AIC for NB and ZINB ------------------- ##
if (file.exists(paste0(distrdir, "/", dataset, exts, "_distribution_fit_summary_data.rds"))) {
  Y <- readRDS(paste0(distrdir, "/", dataset, exts, "_distribution_fit_summary_data.rds"))$GOF_res %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::mutate(tmp = nbinom_standard_aic - zifnbinom_standard_aic) %>%
    dplyr::mutate(asinh_nb_minus_zinb_aic = asinh(tmp)) %>%
    dplyr::mutate(ncr = paste0(ncells, ".", repl)) %>%
    dplyr::select(gene, ncr, asinh_nb_minus_zinb_aic) %>%
    tidyr::spread(ncr, asinh_nb_minus_zinb_aic)
  colnames(Y)[colnames(Y) != "gene"] <- paste0("asinh_nb_minus_zinb_aic.", 
                                               colnames(Y)[colnames(Y) != "gene"])
  truth <- dplyr::full_join(truth, Y, by = "gene")
}

## ---------------------- Define relative truths ---------------------------- ##
## Define "truth" for each method as the genes that are differentially 
## expressed with the largest sample size
(maxn <- max(as.numeric(as.character(get_nsamples(colnames(padj(cobra)))))))
tmp <- padj(cobra)[, get_nsamples(colnames(padj(cobra))) == maxn, drop = FALSE]
colnames(tmp) <- paste0(get_method(colnames(tmp)), ".truth")
tmp <- (tmp <= 0.05)
mode(tmp) <- "numeric"
tmp <- tmp[match(truth$gene, rownames(tmp)), , drop = FALSE]
rownames(tmp) <- truth$gene
## At this point, truth (and thus also tmp) contains all genes that are tested 
## for at least one subset. For the unfiltered data, all genes are tested for 
## the largest subset. However, for the TPM-filtered data, this is not 
## necessarily the case. Despite this, we still consider all genes included in 
## the truth table as "eligible" for relative performance comparisons. Thus, we 
## set NA p-values to 0, since these represent genes that were indeed tested (in
## at least one subset), but filtered out internally
tmp[is.na(tmp)] <- 0

truth <- dplyr::full_join(truth, 
                          data.frame(tmp, stringsAsFactors = FALSE) %>% 
                            dplyr::mutate(gene = rownames(tmp)), 
                          by = "gene")
rownames(truth) <- truth$gene

cobra <- COBRAData(truth = truth, object_to_extend = cobra)

## ------------ Summarize number of calls for each method ------------------- ##
## Make data frame with number of significant, non-significant and NA calls for
## each method. This is calculated based on an adjusted p-value threshold of 0.05.
cobraperf <- calculate_performance(cobra, aspects = "fpr", 
                                   binary_truth = paste0(demethods[1],
                                                         exts, ".truth"), thrs = 0.05)
sign_0.05 <- fpr(cobraperf) %>% dplyr::select(method, NBR, TOT_CALLED) %>%
  dplyr::rename(nbr_sign_adjp0.05 = NBR) %>% dplyr::rename(nbr_called = TOT_CALLED) %>%
  dplyr::mutate(nbr_nonsign_adjp0.05 = nbr_called - nbr_sign_adjp0.05) %>%
  tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.")
tested <- 
  data.frame(nbr_tested = 
               colSums(truth(cobra)[, grep("tested", colnames(truth(cobra)))], 
                       na.rm = TRUE))
tested$dataset <- rownames(tested)
tested <- tested %>% 
  tidyr::separate(dataset, into = c("type", "ncells", "repl"), sep = "\\.")
nbr_called <- dplyr::full_join(sign_0.05, tested) %>% dplyr::select(-type)
nbr_called$dataset <- dataset
nbr_called$filt <- filt
nbr_called$nbr_NA <- nbr_called$nbr_tested - nbr_called$nbr_called

## --------------------- Save output files ---------------------------------- ##
saveRDS(timings, file = paste0(outdir, "/", dataset, exts, "_timings.rds"))
saveRDS(runstatus, file = paste0(outdir, "/", dataset, exts, "_runstatus.rds"))
saveRDS(nbr_called, file = paste0(outdir, "/", dataset, exts, "_nbr_called.rds"))
saveRDS(cobra, file = paste0(outdir, "/", dataset, exts, "_cobra.rds"))

sessionInfo()
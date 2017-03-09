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
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")
source("/home/Shared/data/seq/conquer/comparison/scripts/calculate_gene_characteristics.R")

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

demethods <- strsplit(demethods, ",")[[1]]

print(demethods)
print(dataset)
print(config_file)
print(filt)

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

## Create iCOBRA object from the result files for the different methods
(resfiles <- paste0("/home/Shared/data/seq/conquer/comparison/results/",
                    dataset, "_", demethods, exts, ".rds"))
file.exists(resfiles)
cobra <- NULL
timings <- list()
for (rfn in resfiles) {
  rf <- readRDS(rfn)
  for (nm in names(rf)) {
    print(names(rf[[nm]]))
    timings[[nm]] <- rf[[nm]]$timing
    df <- rf[[nm]]$df
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
    if ("score" %in% colnames(df)) {
      cobra <- COBRAData(score = setNames(data.frame(mt = df$score,
                                                     row.names = rownames(df)), nm),
                         object_to_extend = cobra)
    }
  }
}

cobra <- calculate_adjp(cobra)

## Add gene characteristics to the COBRA object
config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)

subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition

sizes <- names(keep_samples)
truth <- list()
#tested <- list()
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, filt = filt)
    chars <- calculate_gene_characteristics(L)
    df2 <- chars$characs
    colnames(df2)[colnames(df2) != "gene"] <- paste0(colnames(df2)[colnames(df2) != "gene"],
                                                     ".", sz, ".", i)
    truth[[paste0(sz, ".", i)]] <- df2
  }
}
truth <- Reduce(function(...) merge(..., by = "gene", all = TRUE), truth)

## Define "truth" for each method as the genes that are differentially 
## expressed with the largest sample size
tmp <- padj(cobra)[, get_nsamples(colnames(padj(cobra))) == 
                     max(as.numeric(as.character(get_nsamples(colnames(padj(cobra))))))]
colnames(tmp) <- paste0(get_method(colnames(tmp)), ".truth")
tmp <- (tmp <= 0.05)
mode(tmp) <- "numeric"
tmp <- tmp[match(truth$gene, rownames(tmp)), ]
tmp[is.na(tmp)] <- 0

truth <- merge(truth, tmp, by.x = "gene", by.y = 0, all = TRUE)
rownames(truth) <- truth$gene

cobra <- COBRAData(truth = truth, object_to_extend = cobra)

## Make data frame with number of significant, non-significant and NA calls for each method
cobraperf <- calculate_performance(cobra, aspects = "fpr", 
                                   binary_truth = paste0(demethods[1], ".truth"), thrs = 0.05)
sign_0.05 <- fpr(cobraperf) %>% dplyr::select(method, NBR, TOT_CALLED) %>%
  dplyr::rename(nbr_sign0.05 = NBR) %>% dplyr::rename(nbr_called = TOT_CALLED) %>%
  dplyr::mutate(nbr_nonsign0.05 = nbr_called - nbr_sign0.05) %>%
  tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.")
tested <- data.frame(nbr_tested = colSums(truth(cobra)[, grep("tested", colnames(truth(cobra)))], na.rm = TRUE))
tested$dataset <- rownames(tested)
tested <- tested %>% tidyr::separate(dataset, into = c("type", "ncells", "repl"), sep = "\\.")
nbr_called <- inner_join(sign_0.05, tested) %>% dplyr::select(-type)
nbr_called$dataset <- dataset
nbr_called$filt <- filt
nbr_called$nbr_NA <- nbr_called$nbr_tested - nbr_called$nbr_called

saveRDS(cobra, file = paste0("figures/cobra_data/", dataset, exts, "_cobra.rds"))
saveRDS(timings, file = paste0("figures/cobra_data/", dataset, exts, "_timings.rds"))
saveRDS(nbr_called, file = paste0("figures/cobra_data/", dataset, exts, "_nbr_called.rds"))

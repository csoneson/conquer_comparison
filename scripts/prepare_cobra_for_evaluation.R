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
for (rf in resfiles) {
  rf <- readRDS(rf)
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
tested <- list()
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, filt = filt)
    avecount <- data.frame(avecount = apply(L$count, 1, mean), gene = rownames(L$count))
    avetpm <- data.frame(avetpm = apply(L$tpm, 1, mean), gene = rownames(L$tpm))
    fraczero <- data.frame(fraczero = apply(L$count, 1, function(x) mean(x == 0)),
                           fraczero1 = apply(L$count[, L$condt == levels(factor(L$condt))[1]], 
                                             1, function(x) mean(x == 0)),
                           fraczero2 = apply(L$count[, L$condt == levels(factor(L$condt))[2]], 
                                             1, function(x) mean(x == 0)), gene = rownames(L$count))
    fraczero$fraczerodiff <- abs(fraczero$fraczero1 - fraczero$fraczero2)
    vartpm <- data.frame(vartpm = apply(L$tpm, 1, var), gene = rownames(L$tpm))
    df2 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
                  list(vartpm, fraczero, avecount, avetpm))
    colnames(df2)[colnames(df2) != "gene"] <- paste0(colnames(df2)[colnames(df2) != "gene"],
                                                     ".", sz, ".", i)
    truth[[paste0(sz, ".", i)]] <- df2
    
    df3 <- data.frame(gene = df2$gene, tested = TRUE)
    colnames(df3)[colnames(df3) != "gene"] <- paste0(colnames(df3)[colnames(df3) != "gene"],
                                                     ".", sz, ".", i)
    tested[[paste0(sz, ".", i)]] <- df3
  }
}
truth <- Reduce(function(...) merge(..., by = "gene", all = TRUE), truth)
tested <- Reduce(function(...) merge(..., by = "gene", all = TRUE), tested)

padjm <- reshape2::melt(as.matrix(padj(cobra)), value.name = "padj",
                        varnames = c("gene", "method")) %>%
  tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.")
truthm <- reshape2::melt(truth) %>% 
  tidyr::separate(variable, into = c("measurement", "ncells", "repl"), sep = "\\.")
testedm <- reshape2::melt(tested, id.vars = "gene") %>%
  tidyr::separate(variable, into = c("measurement", "ncells", "repl"), sep = "\\.") %>%
  dplyr::rename(tested = value) %>%
  dplyr::select(-measurement)
summary_data <- list(all_data = dplyr::inner_join(padjm, inner_join(truthm, testedm)) %>%
                       dplyr::mutate(dataset = dataset, filt = filt) %>%
                       mutate(value = ifelse(measurement %in% c("vartpm", "avecount"), 
                                             log2(value), value)) %>% 
                       mutate(measurement = ifelse(measurement %in% c("vartpm", "avecount"), 
                                                   paste0("log2_", measurement), measurement)) %>%
                       dplyr::mutate(tested = ifelse(tested == TRUE, TRUE, FALSE)))

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

saveRDS(cobra, file = paste0("figures/cobra_data/", dataset, exts, ".rds"))
saveRDS(timings, file = paste0("figures/cobra_data/", dataset, exts, "_timings.rds"))
saveRDS(summary_data, file = paste0("figures/cobra_data/", dataset, exts, "_summary_data.rds"))
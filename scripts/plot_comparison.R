args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

demethods <- strsplit(demethods, ",")[[1]]

print(demethods)
print(dataset)
print(config_file)
print(filt)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_functions.R")
suppressPackageStartupMessages(library(iCOBRA))

cols <- c("#488d00", "#6400a6", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#a3c6ff", "#e6a900", "#a996ff", "#401600",
          "#ff6d9b", "#017671", "cyan", "red", "blue", "orange",
          "#B17BA6")
if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                        "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                        "voomlimma", "Wilcoxon", "BPSC", "MASTcounts", "MASTcountsDetRate", 
                        "MASTtpm", "zingeRauto"), exts)

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

# pval(cobra)[is.na(pval(cobra))] <- 1
# padj(cobra)[is.na(padj(cobra))] <- 1

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

pdf(paste0("figures/comparison/", dataset, exts, ".pdf"), width = 14, height = 9)
summary_data <- plot_results_relativetruth(cobra, colvec = cols, summary_data = list())
summary_data <- plot_results_relativetruth_all(cobra, colvec = cols, summary_data = summary_data)
summary_data <- plot_results_characterization(cobra, colvec = cols, summary_data = summary_data)
summary_data <- plot_results(cobra, colvec = cols, summary_data = summary_data)
summary_data <- plot_ks(cobra, colvec = cols, summary_data = summary_data)
summary_data <- plot_truefpr(cobra, colvec = cols, summary_data = summary_data)
summary_data <- plot_timing(timings, colvec = cols, summary_data = summary_data)
dev.off()

## Save summary data
summary_data$stats_charac$dataset <- dataset
summary_data$stats_charac$filt <- filt
summary_data$fracpbelow0.05$dataset <- dataset
summary_data$fracpbelow0.05$filt <- filt
summary_data$timing$dataset <- dataset
summary_data$timing$filt <- filt

saveRDS(summary_data, file = paste0("figures/summary_data/", dataset, exts, "_summary_data.rds"))

sessionInfo()

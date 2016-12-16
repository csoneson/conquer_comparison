args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

demethods <- strsplit(demethods, ",")[[1]]

print(demethods)
print(dataset)
print(config_file)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_functions.R")
suppressPackageStartupMessages(library(iCOBRA))

cols <- c("#488d00", "#6400a6", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#a3c6ff", "#e6a900", "#a996ff", "#401600",
          "#ff6d9b", "#017671", "cyan", "red", "blue", "orange", "black")
names(cols) <- c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                 "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                 "voomlimma", "Wilcoxon", "BPSC", "MASTcounts", "MASTcountsDetRate", "MASTtpm")

## Create iCOBRA object
(resfiles <- paste0("/home/Shared/data/seq/conquer/comparison/results/",
                    dataset, "_", demethods, ".rds"))
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

pval(cobra)[is.na(pval(cobra))] <- 1
padj(cobra)[is.na(padj(cobra))] <- 1

pdf(paste0("figures/comparison/", dataset, ".pdf"), width = 14)
plot_results_relativetruth(cobra, colvec = cols)
plot_results_relativetruth_all(cobra, colvec = cols)
plot_results_characterization(cobra, config = config_file, colvec = cols)
plot_results(cobra, colvec = cols)
plot_timing(timings, colvec = cols)
dev.off()

sessionInfo()

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

demethods <- strsplit(demethods, ",")[[1]]

print(demethods)
print(dataset)
print(filt)

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_functions.R")
suppressPackageStartupMessages(library(iCOBRA))

cobras <- list()

cols <- c("#488d00", "#6400a6", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#a3c6ff", "#e6a900", "#a996ff", "#401600",
          "#ff6d9b", "#017671", "cyan", "red", "blue", "orange",
          "#777777", "#7BAFDE", "#F6C141", "#90C987", "#1965B0",
          "#882E72", "#F7EE55")
if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}
names(cols) <- paste0(c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                        "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                        "voomlimma", "Wilcoxon", "BPSC", "MASTcounts", "MASTcountsDetRate", 
                        "MASTtpm", "zingeRauto", "Seurat", "DESeq2census", "edgeRLRTcensus",
                        "DESeq2nofilt", "Seuratnofilt", "NODESnofilt"), exts)

for (tp in c("", "mock")) {
  ## Create iCOBRA object from the result files for the different methods
  (resfiles <- paste0("/home/Shared/data/seq/conquer/comparison/results/",
                      dataset, tp, "_", demethods, exts, ".rds"))
  file.exists(resfiles)
  cobra <- NULL
  timings <- list()
  for (rf in resfiles) {
    print(rf)
    rf <- readRDS(rf)
    for (nm in names(rf)) {
      print(paste0(nm, ":", paste(names(rf[[nm]]), collapse = ", ")))
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
  
  cobras[[paste0("tp_", tp)]] <- cobra
}

pdf(paste0("figures/orig_vs_mock/", dataset, exts, "_orig_vs_mock.pdf"), width = 14, height = 9)
summary_data <- plot_compare_orig_mock(cobras, colvec = cols, summary_data = list())
dev.off()

summary_data$spearman$dataset <- dataset
summary_data$spearman$filt <- filt
summary_data$jaccard_adjp0.05$dataset <- dataset
summary_data$jaccard_adjp0.05$filt <- filt
summary_data$jaccard_adjp0.1$dataset <- dataset
summary_data$jaccard_adjp0.1$filt <- filt
summary_data$jaccard_adjp0.2$dataset <- dataset
summary_data$jaccard_adjp0.2$filt <- filt

saveRDS(summary_data, file = paste0("figures/summary_data/", dataset, exts,
                                    "_summary_data_orig_vs_mock.rds"))

sessionInfo()
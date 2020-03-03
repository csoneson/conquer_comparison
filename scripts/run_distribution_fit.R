args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(powsim))
source("scripts/prepare_mae.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

print(dataset)
print(config_file)
print(filt)
print(figdir)

config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)

subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition
sizes <- names(keep_samples)

(sz <- max(as.numeric(as.character(sizes))))
i <- 1

L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i, 
                imposed_condition = imposed_condition, filt = filt,
                impute = config$impute)

distfit <- powsim::evaluateDist(cnts = round(L$count), RNAseq = "singlecell", 
                                ncores = 1, nsims = 1, frac.genes = 1, 
                                min.meancount = 0.1, min.libsize = 1000)

pdf(paste0(figdir, "/", dataset, exts, "_distribution_fit_summary_data.pdf"), 
    width = 10, height = 10)
powsim::plotEvalDist(evalDist = distfit, annot = FALSE)
dev.off()

distfit$GOF_res$dataset <- dataset
distfit$GOF_res$filt <- filt
distfit$GOF_res$ncells <- sz
distfit$GOF_res$repl <- i

saveRDS(distfit, file = paste0(figdir, "/", dataset, exts, "_distribution_fit_summary_data.rds"))

sessionInfo()
date()

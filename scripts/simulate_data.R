args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(config_file)
print(nDE)
print(seed)

set.seed(seed)

source("/home/Shared/data/seq/conquer/comparison/software/zingeR/R/methods.R")
source("/home/Shared/data/seq/conquer/comparison/software/zingeR/R/simulation.R")

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(gamlss.tr))
suppressPackageStartupMessages(library(rjson))
gamlss.tr::gen.trun(par = 0, family = "NBI", name = "ZeroTruncated",
                    type = "left", varying = FALSE)

config <- fromJSON(file = config_file)

print(config)

mae <- readRDS(config$mae)
pdata <- Biobase::pData(mae)

groupid <- config$groupid
keepgroups <- config$keepgroups
if (length(groupid) > 1) {
  pdata[, paste(groupid, collapse = ".")] <- 
    as.character(interaction(as.data.frame(pdata[, groupid])))
  groupid <- paste(groupid, collapse = ".")
}

counts <- assays(experiments(mae)[["gene"]])[["count_lstpm"]]
stopifnot(all(colnames(counts) == rownames(pdata)))

keep <- which(pdata[, groupid] %in% keepgroups)
length(keep)
counts <- round(counts[, keep])
counts <- counts[rowSums(counts > 0) > 1, ]
group <- as.character(pdata[keep, groupid])

params <- getDatasetZTNB(counts = counts, design = model.matrix(~group))
DEind <- sample(1:nrow(counts), nDE, replace = FALSE)
fcSim <- (2 + rexp(length(DEind), rate = 1/2)) #fold changes

sim <- NBsimSingleCell(dataset = counts, group = factor(group), nTags = nrow(counts), ind = DEind, 
                       params = params, foldDiff = fcSim, lib.size = colSums(counts), verbose = TRUE)

sim_counts <- sim$counts
colnames(sim_counts) <- paste0("s", 1:ncol(sim_counts))
sim_group <- sim$group
sim_foldchange <- sim$foldDiff
sim_indDE <- sim$indDE
sim_propzero <- sim$propZeroGene
sim_dispersion <- sim$Dispersion[, 1]


## Get TPMs from counts. Ideally, we would like to do:
## newMat <- sim_counts/rowMeans(avetxlength)
## sim_TPM <- t(t(newMat) / colSums(newMat)) * 1e6
## However, we don't know which row in the simulated data is obtained from 
## which row in the original data, so we approximate the TPMs by just 
## scaling the (length-scaled TPM) counts to sum to 1e6
sim_TPM <- t(t(sim_counts)/colSums(sim_counts)) * 1e6

generse <- SummarizedExperiment(assays = list(TPM = sim_TPM,
                                              count_lstpm = sim_counts))

## Generate MultiAssayExperiment
mae <- MultiAssayExperiment(experiments = list(gene = generse),
                            pData = data.frame(group = as.character(sim_group),
                                               row.names = colnames(sim_counts),
                                               stringsAsFactors = FALSE))

truth <- rep(0, nrow(sim_counts))
truth[DEind] <- 1
truth <- data.frame(status = truth, dispersion = sim_dispersion, 
                    propzero = sim_propzero, foldchange = sim_foldchange, 
                    stringsAsFactors = FALSE)
rownames(truth) <- rownames(sim_counts)

saveRDS(mae, file = paste0("/home/Shared/data/seq/conquer/comparison/data/", id, "sim", seed, ".rds"))
saveRDS(truth, file = paste0("/home/Shared/data/seq/conquer/comparison/data/", id, "sim", seed, "_truth.rds"))

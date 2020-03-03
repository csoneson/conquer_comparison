args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(dataset)
print(config_file)
print(pDE)
print(seed)

set.seed(seed)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(powsim))

source("scripts/powsim_modified_functions.R")

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

param <- estimateNBParamMod(countData = counts, 
                            cData = data.frame(group = group, sample = colnames(counts), 
                                               row.names = colnames(counts),
                                               stringsAsFactors = FALSE), 
                            design = ~group, RNAseq = "singlecell", 
                            estFramework = "edgeR", sigma = 1.96)

## Check that we keep all genes for the simulation (otherwise the ID matching
## below won't be correct)
length(param$means)
nrow(counts)

## Function to generate log-fold changes
lfc <- function(x) sample(c(-1, 1), size = x, replace = TRUE) * rgamma(x, 4, 2)

desetup <- DESetup(ngenes = nrow(counts), nsims = 1, p.DE = pDE, LFC = lfc, sim.seed = seed)
simsetup <- SimSetup(desetup = desetup, params = param, size.factors = "given")
simsetup$ncores <- 5
simsetup$DEid <- simsetup$DEid[[1]]
simsetup$lfc <- simsetup$lfc[[1]]
simsetup$sim.seed <- simsetup$sim.seed[[1]]

## Simulate
sims <- powsim:::.simRNAseq(simOptions = simsetup, 
                            n1 = sum(group == levels(factor(group))[1]), 
                            n2 = sum(group == levels(factor(group))[2]))

## Figure out to which gene each simulated gene corresponds (replicate some of
## the sampling code used in the simulation function)
set.seed(simsetup$sim.seed)
index <- sample(1:length(simsetup$means), size = simsetup$ngenes, replace = TRUE)
cor(apply(counts[index, ], 1, mean), apply(sims$counts, 1, mean), method = "spearman")
cor(simsetup$means[index], apply(sims$counts, 1, mean), method = "spearman")
idconv <- data.frame(id = rownames(sims$counts), gene = rownames(counts)[index], 
                     stringsAsFactors = FALSE)

## Generate truth
truth <- data.frame(id = rownames(sims$counts),
                    lfc = sims$simOptions$lfc, stringsAsFactors = FALSE) %>%
  dplyr::mutate(status = as.numeric(lfc != 0))
rownames(truth) <- truth$id

## Calculate TPMs from the simulated counts
avetxlength <- assays(experiments(mae)[["gene"]])[["avetxlength"]]
avetxlength <- avetxlength[match(idconv$gene[match(rownames(sims$counts), 
                                                   idconv$id)], rownames(avetxlength)), ]
newMat <- sims$counts/rowMeans(avetxlength)
sim_TPM <- t(t(newMat) / colSums(newMat)) * 1e6

## Generate synthetic avetxlength matrix
sim_avetxlength <- matrix(rep(rowMeans(avetxlength), each = ncol(sims$counts)), 
                          ncol = ncol(sims$counts), byrow = TRUE)
rownames(sim_avetxlength) <- rownames(sims$counts)
colnames(sim_avetxlength) <- colnames(sims$counts)

## Generate MultiAssayExperiment
generse <- SummarizedExperiment(assays = list(TPM = sim_TPM,
                                              count_lstpm = sims$counts,
                                              avetxlength = sim_avetxlength))

mae <- MultiAssayExperiment(experiments = list(gene = generse),
                            pData = data.frame(group = as.character(levels(factor(group))[sims$designs/2 + 3/2]),
                                               sample = colnames(sims$counts), 
                                               row.names = colnames(sims$counts)))

saveRDS(mae, file = paste0("data/", dataset, "sim", seed, ".rds"))
saveRDS(truth, file = paste0("data/", dataset, "sim", seed, "_truth.rds"))

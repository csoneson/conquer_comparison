## Read Usoskin data
x <- read.delim("data/Usoskin_External_resources_Table_1.txt", header = FALSE, as.is = TRUE)

sample_annot <- x[1:10, -(1:9)]
sample_annot <- t(sample_annot)
colnames(sample_annot) <- x[1:10, 9]
sample_annot <- data.frame(sample_annot, stringsAsFactors = FALSE)

gene_annot <- x[-(1:11), 1:8]
colnames(gene_annot) <- x[11, 1:8]
gene_annot <- data.frame(gene_annot, stringsAsFactors = FALSE)

rpms <- as.matrix(x[-(1:11), -(1:9)])
mode(rpms) <- "numeric"
rownames(rpms) <- x[-(1:11), 1]
colnames(rpms) <- x[1, -(1:9)]

## Filter out empty wells
keep <- which(sample_annot$Content == "cell")
sample_annot <- sample_annot[keep, ]
rpms <- rpms[, keep]

# Filter on cell type
keep <- which(sample_annot$Level.1 %in% c("NF", "NP", "PEP", "TH"))
sample_annot <- sample_annot[keep, ]
rpms <- rpms[, keep]

counts <- sweep(rpms, 2, as.numeric(as.character(sample_annot$Reads))/1e6, "*")

sample_annot$Picking.sessions.Level.3 <- interaction(sample_annot$Picking.sessions, 
                                                     sample_annot$Level.3)
rownames(sample_annot) <- sample_annot$Sample.ID

## RT-1.NP1 vs RT-1.TH

library(SummarizedExperiment)
library(MultiAssayExperiment)

generse <- SummarizedExperiment(assays = list(TPM = rpms,
                                              count = counts,
                                              count_lstpm = counts))

## Generate MultiAssayExperiment
mae <- MultiAssayExperiment(experiments = list(gene = generse),
                            pData = droplevels(sample_annot))
mae@metadata <- list(organism = "Mus musculus")

saveRDS(mae, file = paste0("data/UsoskinGSE59739.rds"))

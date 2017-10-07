suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(akima))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(lazyeval))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(pheatmap))
source("scripts/prepare_mae.R")

#' Convenience functions to extract the first, second and third part of each of
#' a vector of strings (parts are separated by .)
#' 
get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

cols <- c("#771155", "#AA4488", "#CC99BB", 
          "#3A0027", "#BCA2B4",
          "#000000", "#4477AA", "#BEBEBE", 
          "#117777", "#44AAAA", "#77CCCC",
          "#117744", "#44AA77", "#88CCAA", "#004120",
          "#777711", "#AAAA44", "#C1C185", 
          "#774411", "#AA7744", "#827262", "#DBCDBE",
          "#771122", "#AA4455", "#DD7788", "#FFFF00",
          "#DDDD77", "#DDAA77", "#86861B",
          "#114477", "#77AADD",
          "#FFA500", "#7CFC00",
          "#00BFFF", "#00688B",
          "#FF99ED", "#3B5D93")
names(cols) <- c("edgeRLRT", "edgeRLRTdeconv", "edgeRLRTcensus", 
                 "edgeRQLF", "edgeRLRTrobust",
                 "ROTSvoom", "ROTScpm", "ROTStpm",
                 "DESeq2", "DESeq2nofilt", "DESeq2census",
                 "MASTtpm", "MASTcpm", "MASTtpmDetRate", "MASTcpmDetRate",
                 "monocle", "monoclecensus", "monoclecount", 
                 "SeuratBimod", "SeuratBimodnofilt", "SeuratTobit", "SeuratBimodIsExpr2", 
                 "SAMseq", "Wilcoxon", "NODES", "NODESnofilt",
                 "SCDE", "BPSC", "D3E", 
                 "voomlimma", "limmatrend",
                 "metagenomeSeq", "scDD",
                 "zingeRedgeR", "zingeRedgeRnofilt",
                 "ttest", "DEsingle")

## Write color definitions to file
s <- sapply(1:length(cols), function(x) {
  paste0("\\", "definecolor{", names(cols)[x], "}{HTML}{", gsub("#", "", cols[x]), "}")
})
write.table(s, file = "color_definitions.txt", quote = FALSE, 
            sep = "\n", row.names = FALSE, col.names = FALSE)

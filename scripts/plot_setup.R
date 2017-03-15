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
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")

#' Convenience functions to extract the first, second and third part of each of
#' a vector of strings (parts are separated by .)
#' 
get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

cols <- c("#488d00", "black", "#8bff58", "#ff5cd5", "#9CC0AD",
          "#ab0022", "#7c9e58", "#e6a900", "#ff516e", "#364922",
          "#db0a4c", "#017671", "cyan", "#6d6de8", "blue",
          "#a5a5e5", "#777777", "#6400a6", "#f2c6cf", "#afeda6", 
          "#af5d6d", "#bf8bb2", "#F7EE55", "gray", 
          "pink", "#7BAFDE", "#42425b")
names(cols) <- c("edgeRLRT", "zingeR", "SAMseq", "edgeRQLF", "NODES",
                 "DESeq2", "edgeRLRTdeconv", "SCDE", "monocle", "edgeRLRTrobust", 
                 "voomlimma", "Wilcoxon", "BPSC", "MASTcpm", "MASTcpmDetRate", 
                 "MASTtpm", "zingeRauto", "Seurat", "DESeq2census", "edgeRLRTcensus",
                 "DESeq2nofilt", "Seuratnofilt", "NODESnofilt", "zingeRautonofilt",
                 "monoclecensus", "D3E", "MASTtpmDetRate")


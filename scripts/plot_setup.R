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
          "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC",
          "#117744", "#44AA77", "#88CCAA", "#004120",
          "#777711", "#AAAA44", 
          "#774411", "#AA7744",  
          "#771122", "#AA4455", "#DD7788", "yellow",
          "#DDDD77", "#DDAA77", "#86861B",
          "black", "grey")
names(cols) <- c("edgeRLRT", "edgeRLRTdeconv", "edgeRLRTcensus", 
                 "edgeRQLF", "edgeRLRTrobust",
                 "zingeR", "zingeRauto", "zingeRautonofilt",
                 "DESeq2", "DESeq2nofilt", "DESeq2census",
                 "MASTtpm", "MASTcpm", "MASTtpmDetRate", "MASTcpmDetRate",
                 "monocle", "monoclecensus", 
                 "Seurat", "Seuratnofilt", 
                 "SAMseq", "Wilcoxon", "NODES", "NODESnofilt",
                 "SCDE", "BPSC", "D3E", 
                 "voomlimma", "limmatrend")



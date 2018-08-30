suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(forcats))
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

## Define plot characters for the different sample sizes
pch_ncells <- c(16, 17, 15, 3, 7, 8, 4, 6, 9, 10, 11, 12, 13, 14, 1, 2, 5, 18)
names(pch_ncells) <- as.character(c(6, 12, 22, 24, 26, 36, 44, 47, 48, 50, 58, 
                                    90, 95, 96, 120, 200, 240, 400))

## Define colors for the different methods
cols <- c("#771155", "#AA4488", "#CC99BB", 
          "#3A0027", "#BCA2B4", "#FF02AC",
          "#000000", "#4477AA", "#BEBEBE", 
          "#117777", "#44AAAA", "#77CCCC", "#77D6B2", "#21FFFF", 
          "#117744", "#44AA77", "#88CCAA", "#004120",
          "#777711", "#AAAA44", "#C1C185", 
          "#774411", "#AA7744", "#827262", "#DBCDBE",
          "#771122", "#AA4455", "#DD7788", "#FFFF00",
          "#DDDD77", "#DDAA77", "#86861B",
          "#114477", "#77AADD", "#0280FF", 
          "#FFA500", "#7CFC00",
          "#00BFFF", "#00688B",
          "#FF99ED", "#3B5D93",
          "#4286F4", "#98BAF2",
          "#18386D")
names(cols) <- c("edgeRLRT", "edgeRLRTdeconv", "edgeRLRTcensus", 
                 "edgeRQLF", "edgeRLRTrobust", "edgeRQLFDetRate",
                 "ROTSvoom", "ROTScpm", "ROTStpm",
                 "DESeq2", "DESeq2nofilt", "DESeq2census", "DESeq2betapFALSE", "DESeq2LRT", 
                 "MASTtpm", "MASTcpm", "MASTtpmDetRate", "MASTcpmDetRate",
                 "monocle", "monoclecensus", "monoclecount", 
                 "SeuratBimod", "SeuratBimodnofilt", "SeuratTobit", "SeuratBimodIsExpr2", 
                 "SAMseq", "Wilcoxon", "NODES", "NODESnofilt",
                 "SCDE", "BPSC", "D3E", 
                 "voomlimma", "limmatrend", "limmatrendDetRate",
                 "metagenomeSeq", "scDD",
                 "zingeRedgeR", "zingeRedgeRnofilt",
                 "ttest", "DEsingle",
                 "zinbwaveedgeR", "zinbwaveDESeq2",
                 "logregLRT")

## Write color definitions to file
s <- sapply(1:length(cols), function(x) {
  paste0("\\", "definecolor{", names(cols)[x], "}{HTML}{", gsub("#", "", cols[x]), "}")
})
write.table(s, file = "color_definitions.txt", quote = FALSE, 
            sep = "\n", row.names = FALSE, col.names = FALSE)

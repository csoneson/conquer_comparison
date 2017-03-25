args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

print(outdir)
print(figdir)

## --------------------- Make figures for publication ----------------------- ##

## fraction NA
xorig <- readRDS(paste0(figdir, "/fracNA/summary_fracNA.rds"))
xfilt <- readRDS(paste0(figdir, "/fracNA/summary_fracNA_TPM_1_25p.rds"))

## True FPR

## True FDR

## True TPR

## Timing

sessionInfo()
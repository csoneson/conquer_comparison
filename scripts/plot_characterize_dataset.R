args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

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
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(lazyeval))
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

pdf(paste0("figures/dataset_characteristics/", dataset, exts, ".pdf"), width = 14, height = 9)

print(dataset)
print(config_file)
print(filt)

config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)

subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition

sizes <- names(keep_samples)
char_gene <- list()
char_cells <- list()
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, filt = filt)
    
    ## Count distributions
    print(reshape2::melt(L$count) %>% ggplot(aes(x = value, group = Var2)) + 
            geom_density() + scale_x_log10() + theme_bw() + 
            xlab("Count") + ggtitle(paste0("Count distribution per cell, ", 
                                           sz, " cells per group, repl ", i)))
    
    ## Gene characteristics
    avecount <- data.frame(avecount = apply(L$count, 1, mean), gene = rownames(L$count))
    avetpm <- data.frame(avetpm = apply(L$tpm, 1, mean), gene = rownames(L$tpm))
    fraczero <- data.frame(fraczero = apply(L$count, 1, function(x) mean(x == 0)),
                           fraczero1 = apply(L$count[, L$condt == levels(factor(L$condt))[1]], 
                                             1, function(x) mean(x == 0)),
                           fraczero2 = apply(L$count[, L$condt == levels(factor(L$condt))[2]], 
                                             1, function(x) mean(x == 0)), gene = rownames(L$count))
    fraczero$fraczerodiff <- abs(fraczero$fraczero1 - fraczero$fraczero2)
    vartpm <- data.frame(vartpm = apply(L$tpm, 1, var), gene = rownames(L$tpm))
    df2 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
                  list(vartpm, fraczero, avecount, avetpm))
    colnames(df2)[colnames(df2) != "gene"] <- paste0(colnames(df2)[colnames(df2) != "gene"],
                                                     ".", sz, ".", i)
    char_gene[[paste0(sz, ".", i)]] <- df2
    
    ## Cell characteristics
    libsize <- data.frame(libsize = colSums(L$count), cell = colnames(L$count))
    fraczerocell <- data.frame(fraczero = colMeans(L$count == 0), cell = colnames(L$count))
    df3 <- Reduce(function(...) merge(..., by = "cell", all = TRUE),
                  list(libsize, fraczerocell))
    colnames(df3)[colnames(df3) != "cell"] <- paste0(colnames(df3)[colnames(df3) != "cell"],
                                                     ".", sz, ".", i)
    char_cells[[paste0(sz, ".", i)]] <- df3
  }
}
char_gene <- Reduce(function(...) merge(..., by = "gene", all = TRUE), char_gene)
char_gene_m <- reshape2::melt(char_gene) %>% 
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group"))) 

char_cells <- Reduce(function(...) merge(..., by = "cell", all = TRUE), char_cells)
char_cells_m <- reshape2::melt(char_cells) %>% 
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group"))) 

for (tp in c("vartpm", "avecount", "avetpm")) {
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + scale_x_log10() + facet_wrap(~ncells) + 
          theme_bw() + xlab(ifelse(tp == "vartpm", "Variance of TPM values per gene", 
                                   ifelse(tp == "avecount", "Average count per gene",
                                          "Average TPM per gene"))))
}
for (tp in c("fraczero", "fraczerodiff")) {
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + facet_wrap(~ncells) + 
          theme_bw() + xlab(ifelse(tp == "fraczero", "Fraction zeros per gene", 
                                   "Difference (between conditions) of zero fraction per gene")))
}
for (tp in c("libsize", "fraczero")) {
  print(char_cells_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + facet_wrap(~ncells) + 
          theme_bw() + xlab(ifelse(tp == "fraczero", "Fraction zeros per cell", 
                                   "Library size per cell")))
}

dev.off()

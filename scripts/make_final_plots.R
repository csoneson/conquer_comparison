args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

print(figdir)  ## Directory where individual plots are stored
print(outdir)  ## Directory where final plots will be saved

## --------------------- Make figures for publication ----------------------- ##

## t-SNEs
X <- list()
for (ds in c("GSE74596", "GSE45719", "GSE48968-GPL13112", "EMTAB2805", 
             "UsoskinGSE59739", "GSE63818-GPL16791", "GSE60749-GPL13112", "EGEUV1")) {
  X[[ds]] <- readRDS(paste0(figdir, "/dataset_characteristics/", ds,
                            "_dataset_characteristics_plots.rds"))
}
pdf(paste0(outdir, "/tsne_final.pdf"), width = 9, height = 22)
print(plot_grid(plotlist = lapply(X, function(x) x$tsne + guides(color = guide_legend(nrow = 2))), 
                ncol = 2, labels = LETTERS[1:8], align = "h"))
dev.off()

## fraction NA
xorig <- readRDS(paste0(figdir, "/multi_dataset/fracNA/summary_fracNA_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/fracNA/summary_fracNA_TPM_1_25p_plots.rds"))
p <- plot_grid(plot_grid(xorig$fracna_comb + theme(legend.position = "none") + 
                           ggtitle("Without filtering"), 
                         xfilt$fracna_comb + theme(legend.position = "none") + 
                           ggtitle("After filtering"),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$fracna_comb + 
                            theme(legend.position = "bottom") + 
                            guides(colour = FALSE,
                                   shape = 
                                     guide_legend(nrow = 1,
                                                  title = "Number of cells per group",
                                                  title.theme = element_text(size = 12,
                                                                             angle = 0),
                                                  label.theme = element_text(size = 10,
                                                                             angle = 0),
                                                  keywidth = 1, default.unit = "cm"))),
               rel_heights = c(1.7, 0.1), ncol = 1)
pdf(paste0(outdir, "/fracNA_final.pdf"), width = 12, height = 6)
print(p)
dev.off()

## True FPR
xorig <- readRDS(paste0(figdir, "/multi_dataset/truefpr/summary_truefpr_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/truefpr/summary_truefpr_TPM_1_25p_plots.rds"))
p <- plot_grid(plot_grid(xorig$truefpr + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(0,0.8), 
                         xfilt$truefpr + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(0,0.8),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$truefpr + 
                            theme(legend.position = "bottom") + 
                            guides(colour = FALSE,
                                   shape = 
                                     guide_legend(nrow = 1,
                                                  title = "Number of cells per group",
                                                  title.theme = element_text(size = 12,
                                                                             angle = 0),
                                                  label.theme = element_text(size = 10,
                                                                             angle = 0),
                                                  keywidth = 1, default.unit = "cm"))),
               rel_heights = c(1.7, 0.1), ncol = 1)
pdf(paste0(outdir, "/truefpr_final.pdf"), width = 12, height = 6)
print(p)
dev.off()

## True FDR

## True TPR

## Timing

saveRDS(NULL, file = paste0(outdir, "/final_plots.rds"))

sessionInfo()
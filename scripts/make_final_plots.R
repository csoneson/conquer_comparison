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
pdf(paste0(outdir, "/tsne_final.pdf"), width = 13.5, height = 16.5)
print(plot_grid(plotlist = lapply(X, function(x) x$tsne + guides(color = guide_legend(nrow = 2))), 
                ncol = 3, labels = LETTERS[1:8], align = "h"))
dev.off()

## ---------------------- Result characterization --------------------------- ##
x <- readRDS(paste0(figdir, "/multi_dataset/pca/summary_pca_plots.rds"))
p <- x[["GSE74596mock_GSE45719mock_EMTAB2805mock_GSE60749-GPL13112mock_GSE48968-GPL13112mock_UsoskinGSE59739mocktstat_distribution"]] + 
  theme(legend.position = "bottom") + 
  guides(colour = FALSE,
         shape = guide_legend(nrow = 2,
                              title = "",
                              title.theme = element_text(size = 12,
                                                         angle = 0),
                              label.theme = element_text(size = 12,
                                                         angle = 0)))
pdf(paste0(outdir, "/result_characterization_final.pdf"), width = 10, height = 6)
print(p)
dev.off()

x <- readRDS(paste0(figdir, "/multi_dataset/pca/summary_pca_sim_plots.rds"))
p <- x[["GSE45719sim123mock_GSE74596sim123mock_GSE48968-GPL13112sim123mocktstat_distribution"]] + 
  theme(legend.position = "bottom") + 
  guides(colour = FALSE,
         shape = guide_legend(nrow = 2,
                              title = "",
                              title.theme = element_text(size = 12,
                                                         angle = 0),
                              label.theme = element_text(size = 12,
                                                         angle = 0)))
pdf(paste0(outdir, "/result_characterization_final_sim.pdf"), width = 10, height = 6)
print(p)
dev.off()

## ----------------------------- concordance -------------------------------- ##
xorig <- readRDS(paste0(figdir, "/multi_dataset/orig_vs_mock/summary_orig_vs_mock_plots.rds"))
p <- plot_grid(plot_grid(xorig$auc_signal_comb + theme(legend.position = "none"), 
                         xorig$auc_mock_comb + theme(legend.position = "none"),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               xorig$tstat_auc + theme(legend.position = "none"),
               get_legend(xorig$tstat_auc + theme(legend.position = "bottom") + 
                 guides(colour = FALSE,
                        shape = 
                          guide_legend(nrow = 2,
                                       title = "Data set",
                                       title.theme = element_text(size = 12,
                                                                  angle = 0),
                                       label.theme = element_text(size = 12,
                                                                  angle = 0),
                                       keywidth = 1, default.unit = "cm"))),
               rel_heights = c(1.7, 1.7, 0.1), ncol = 1, labels = c("", "C", ""))
pdf(paste0(outdir, "/concordancediff_final.pdf"), width = 12, height = 12)
print(p)
dev.off()

## ---------------------------- fraction NA --------------------------------- ##
## single-cell
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
                                                  label.theme = element_text(size = 12,
                                                                             angle = 0),
                                                  keywidth = 1, default.unit = "cm"))),
               rel_heights = c(1.7, 0.1), ncol = 1)
pdf(paste0(outdir, "/fracNA_final.pdf"), width = 12, height = 6)
print(p)
dev.off()

## Split by data set
xorig <- readRDS(paste0(figdir, "/multi_dataset/fracNA/summary_fracNA_plots.rds"))
p <- xorig$fracna_sep + theme(legend.position = "bottom") + 
  ggtitle("Without filtering") + 
  guides(colour = FALSE,
         shape = guide_legend(nrow = 1,
                        title = "Number of cells per group",
                        title.theme = element_text(size = 12, angle = 0),
                        label.theme = element_text(size = 12, angle = 0),
                        keywidth = 1, default.unit = "cm"))
pdf(paste0(outdir, "/fracNA_final_sepbyds.pdf"), width = 12, height = 6)
print(p)
dev.off()

## bulk
xorig <- readRDS(paste0(figdir, "/multi_dataset/fracNA/summary_fracNA_bulk_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/fracNA/summary_fracNA_TPM_1_25p_bulk_plots.rds"))
p <- plot_grid(plot_grid(xorig$fracna_comb + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(-0.01, 1), 
                         xfilt$fracna_comb + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(-0.01, 1),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$fracna_comb + 
                            theme(legend.position = "bottom") + 
                            guides(colour = FALSE,
                                   shape = 
                                     guide_legend(nrow = 1,
                                                  title = "Number of cells per group",
                                                  title.theme = element_text(size = 12,
                                                                             angle = 0),
                                                  label.theme = element_text(size = 12,
                                                                             angle = 0),
                                                  keywidth = 1, default.unit = "cm"))),
               rel_heights = c(1.7, 0.1), ncol = 1)
pdf(paste0(outdir, "/fracNA_final_bulk.pdf"), width = 12, height = 6)
print(p)
dev.off()

## ------------------------------ True FPR ---------------------------------- ##
## single-cell
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

## simulated
xorig <- readRDS(paste0(figdir, "/multi_dataset/truefpr/summary_truefpr_sim_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/truefpr/summary_truefpr_TPM_1_25p_sim_plots.rds"))
p <- plot_grid(plot_grid(xorig$truefpr + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(0, 0.8), 
                         xfilt$truefpr + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(0, 0.8),
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
pdf(paste0(outdir, "/truefpr_final_sim.pdf"), width = 12, height = 6)
print(p)
dev.off()

## True FDR
xorig <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_sim_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_TPM_1_25p_sim_plots.rds"))
p <- plot_grid(plot_grid(xorig$FDR_all + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(0,1), 
                         xfilt$FDR_all + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(0,1),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$FDR_all + 
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
pdf(paste0(outdir, "/truefdr_final.pdf"), width = 12, height = 6)
print(p)
dev.off()

## True TPR
xorig <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_sim_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_TPM_1_25p_sim_plots.rds"))
p <- plot_grid(plot_grid(xorig$TPR_all + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(0,1), 
                         xfilt$TPR_all + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(0,1),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$TPR_all + 
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
pdf(paste0(outdir, "/truetpr_final.pdf"), width = 12, height = 6)
print(p)
dev.off()

## True AUROC
xorig <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_sim_plots.rds"))
xfilt <- readRDS(paste0(figdir, "/multi_dataset/trueperformance/summary_trueperformance_TPM_1_25p_sim_plots.rds"))
p <- plot_grid(plot_grid(xorig$auroc_all + theme(legend.position = "none") + 
                           ggtitle("Without filtering") + ylim(0,1), 
                         xfilt$auroc_all + theme(legend.position = "none") + 
                           ggtitle("After filtering") + ylim(0,1),
                         labels = c("A", "B"), align = "h", rel_widths = c(1, 1), nrow = 1),
               get_legend(xorig$auroc_all + 
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
pdf(paste0(outdir, "/trueauroc_final.pdf"), width = 12, height = 6)
print(p)
dev.off()


saveRDS(NULL, file = paste0(outdir, "/final_plots.rds"))

sessionInfo()
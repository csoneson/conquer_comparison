
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



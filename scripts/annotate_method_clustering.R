args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(hclustrds)
print(chartxt)

suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(cowplot))

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 5.5)

annot <- read.delim(chartxt, header = TRUE, as.is = TRUE)
y <- readRDS(hclustrds)
stab <- y$stability_scores
subcl <- y$subclusters_average
clust <- y$clustering_average

colors <- list(
  input = c(counts = "#01368C", TPM = "#2F6CCE", CPM = "#93B8F2", Census = "#D6E4F9"),
  modeling = c(parametric = "#C10359", nonparametric = "#E0B1C6"),
  transformation = c(no = "#2EA801", log = "#96E878"),
  NAresults = c(yes = "#4C4001", no = "#8E8D86")
)

ggt <- ggtree(as.phylo(clust))
for (m in seq_len(length(subcl))) {
  tryCatch({
    i <- MRCA(ggt, subcl[[m]])
    if (stab[m] >= 0.1)
      ggt$data[i, "label"] <- round(stab[m], 2)
  }, error = function(e) NULL)
}
ggt <- ggt + geom_label2(aes(subset = !isTip, label = label), size = 2) + 
  geom_tiplab(aes(angle = 90), hjust = 1) + 
  ggplot2::scale_x_reverse() + ggplot2::coord_flip() + 
  theme(plot.margin = unit(c(0, 0, 10, 0), "mm")) + 
  xlim_tree(0.85)

tiporder <- ggt$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
annot <- annot[match(tiporder$label, annot$method), 
               c("input", "modeling", "transformation", "NAresults"), drop = FALSE]

g1 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = input)) + 
  geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
  scale_fill_manual(values = colors$input) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 7, 0, 7), "mm"))
g2 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = modeling)) + 
  geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
  scale_fill_manual(values = colors$modeling) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 7, 0, 7), "mm"))
g3 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = transformation)) + 
  geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
  scale_fill_manual(values = colors$transformation) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 7, 0, 7), "mm"))
g4 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = NAresults)) + 
  geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
  scale_fill_manual(values = colors$NAresults) + 
  scale_x_continuous(expand = c(0, 0), labels = tiporder$label, breaks = 1:nrow(annot)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 7, 0, 7), "mm"))

gds <- plot_grid(get_legend(g1 + 
                              theme(legend.position = "bottom") + 
                              guides(fill = 
                                       guide_legend(ncol = 2,
                                                    title = "Input", title.position = "top",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 10,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 get_legend(g2 + 
                              theme(legend.position = "bottom") + 
                              guides(fill = 
                                       guide_legend(ncol = 1,
                                                    title = "Modeling", title.position = "top",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 10,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 get_legend(g3 + 
                              theme(legend.position = "bottom") + 
                              guides(fill = 
                                       guide_legend(ncol = 1,
                                                    title = "Transformation", title.position = "top",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 10,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 get_legend(g4 + 
                              theme(legend.position = "bottom") + 
                              guides(fill = 
                                       guide_legend(ncol = 1,
                                                    title = "NA values", title.position = "top",
                                                    override.aes = list(size = 1.5),
                                                    title.theme = element_text(size = 12,
                                                                               angle = 0),
                                                    label.theme = element_text(size = 10,
                                                                               angle = 0),
                                                    keywidth = 1, default.unit = "cm"))),
                 nrow = 1, rel_widths = c(1.25, 1, 1, 1))

plot_grid(ggt + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")), 
          g1 + theme(legend.position = "none"), 
          g2 + theme(legend.position = "none"), 
          g3 + theme(legend.position = "none"), 
          g4 + theme(legend.position = "none"),
          gds + theme(plot.margin = unit(c(0, 0, 0, 15), "mm")), 
          rel_heights = c(6, 0.75, 0.75 ,0.75 ,0.75, 2), ncol = 1)

dev.off()

saveRDS(NULL, outrds)
sessionInfo()
date()


args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(hclustrds)
print(chartxt)
print(outrds)

suppressPackageStartupMessages(library(ComplexHeatmap))

x <- read.delim(chartxt, header = TRUE, as.is = TRUE)
y <- readRDS(hclustrds)
y <- y$average_auc_clustering$tree_row

x <- x[match(y$labels, x$method), 
       c("input", "modeling", "transformation", "NAresults"), drop = FALSE]

colors <- list(
  input = c(counts = "#01368C", TPM = "#2F6CCE", CPM = "#93B8F2", Census = "#D6E4F9"),
  modeling = c(parametric = "#C10359", nonparametric = "#E0B1C6"),
  transformation = c(no = "#2EA801", log = "#96E878"),
  NAresults = c(yes = "#4C4001", no = "#8E8D86")
)

ha = HeatmapAnnotation(df = x,
                       text = column_anno_text(y$labels, rot = 90, just = "right"),
                       show_annotation_name = FALSE,
                       col = colors)
ht = Heatmap(matrix(nrow = 0, ncol = nrow(x)), cluster_columns = y, 
             top_annotation = ha, column_title = "",
             column_dend_height = unit(35, "mm"), top_annotation_height = unit(40, "mm"))
pdf(gsub("rds$", "pdf", outrds), width = 10, height = 4.5)
draw(ht, padding = unit(c(2, 2, -35, 2), "mm"))
dev.off()

saveRDS(NULL, outrds)
sessionInfo()
date()

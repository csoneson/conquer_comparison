args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

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
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(lazyeval))
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(scater))
source("scripts/prepare_mae.R")
source("scripts/calculate_gene_characteristics.R")
source("scripts/calculate_cell_characteristics.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

print(dataset)
print(config_file)
print(filt)
print(cell_cycle_file)
print(figdir)

plots <- list()

config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)

subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition

cell_cycle_file <- readRDS(cell_cycle_file)
if (!is.null(metadata(mae)$index)) {
  cell_cycle_ids <- cell_cycle_file[[gsub("\\.cdna.*", "", metadata(mae)$index)]]
} else if (!is.null(metadata(mae)$organism)) {
  cell_cycle_ids <- cell_cycle_file[[gsub(" ", "_", metadata(mae)$organism)]]
} else {
  cell_cycle_ids <- NULL
}

## Make t-SNE plot, indicating the selected groups. Only for unfiltered data,
## since the filtering will anyway take place later
if (filt == "") {
  pdf(paste0(figdir, "/", dataset, exts, "_dataset_characteristics_1.pdf"), width = 7, height = 7)
  groupidcomb <- paste(groupid, collapse = ".")
  sceset <- newSCESet(countData = assays(experiments(mae)[["gene"]])[["count_lstpm"]], 
                      tpmData = assays(experiments(mae)[["gene"]])[["TPM"]],
                      phenoData = new("AnnotatedDataFrame", data = as.data.frame(pData(mae))))
  sceset <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], 
                             return_SCESet = TRUE, draw_plot = FALSE)
  df <- data.frame(cell = rownames(sceset@reducedDimension), 
                   sceset@reducedDimension) %>% 
    dplyr::full_join(data.frame(cell = rownames(pData(sceset)), 
                                group = pData(sceset)[, groupidcomb]),
                     by = "cell")
  plots[["tsne"]] <- 
    ggplot(df, aes(x = X1, y = X2)) + geom_point(size = 2, color = "grey") + 
    geom_point(data = df %>% dplyr::filter(group %in% config$keepgroups), size = 2, aes(color = group)) + 
    scale_color_manual(
      values = structure(
        c("blue", "red", rep("grey", length(unique(pData(sceset)[, groupidcomb])) - 2)), 
        names = c(config$keepgroups, 
                  setdiff(unique(pData(sceset)[, groupidcomb]), config$keepgroups))),
      breaks = config$keepgroups, name = "") + 
    theme_bw() + xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + 
    theme(legend.position = "bottom",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13)) + 
    ggtitle(dataset)
  print(plots[["tsne"]])
  dev.off()
}

pdf(paste0(figdir, "/", dataset, exts, "_dataset_characteristics_2.pdf"), width = 14, height = 9)

sizes <- names(keep_samples)
char_gene <- list()
char_cells <- list()
char_ds <- list()
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i, 
                    imposed_condition = imposed_condition, filt = filt,
                    impute = config$impute)

    ## Gene characteristics
    chars <- calculate_gene_characteristics(L, do.plot = TRUE, 
                                            title.ext = paste0(", ", sz, " cells per group, repl ", i))
    df2 <- chars$characs
    colnames(df2)[colnames(df2) != "gene"] <- paste0(colnames(df2)[colnames(df2) != "gene"],
                                                     ".", sz, ".", i)
    char_gene[[paste0(sz, ".", i)]] <- df2
    
    ## Cell characteristics
    cellchars <- calculate_cell_characteristics(L)
    df3 <- cellchars$characs
    df3$condition <- L$condt[match(df3$cell, names(L$condt))]
    colnames(df3)[colnames(df3) != "cell"] <- paste0(colnames(df3)[colnames(df3) != "cell"],
                                                     ".", sz, ".", i)
    char_cells[[paste0(sz, ".", i)]] <- df3
    
    ## Dataset characteristics
    char_ds[[paste0(sz, ".", i)]] <- c(n_genes = nrow(L$count),
                                       silhouette_avg = mean(df3[, paste0("silhouette.", sz, ".", i)]))
    
    ## Plot expression of cell cycle genes
    if (!is.null(cell_cycle_ids)) {
      tryCatch({
        tpm_cell_cycle <- L$tpm[match(cell_cycle_ids, rownames(L$tpm)), ]
        tpm_cell_cycle <- tpm_cell_cycle[!is.na(rownames(tpm_cell_cycle)), ]
        tpm_cell_cycle_m <- reshape2::melt(tpm_cell_cycle)
        tpm_cell_cycle_m$condition <- L$condt[match(tpm_cell_cycle_m$Var2, names(L$condt))]
        tpm_cell_cycle_m$Var2 <- factor(tpm_cell_cycle_m$Var2, levels = names(sort(L$condt)))
        nr <- nrow(tpm_cell_cycle)
        vargroup <- data.frame(gene = unique(tpm_cell_cycle_m$Var1), group = rep(1:25, ceiling(nr/25))[1:nr])
        tpm_cell_cycle_m$plot_group <- vargroup$group[match(tpm_cell_cycle_m$Var1, vargroup$gene)]
        tpm_cell_cycle_m <- tpm_cell_cycle_m %>% dplyr::group_by(plot_group) %>% 
          dplyr::mutate(plot_color = paste0("p", as.numeric(as.factor(Var1))))
        print(ggplot(tpm_cell_cycle_m, aes(x = Var2, y = value + 1)) + 
                geom_line(aes(group = Var1, color = plot_color)) + geom_point(aes(shape = condition)) + 
                guides(color = FALSE) + scale_y_log10() + 
                ggtitle(paste0("Cell cycle-associated genes, ", sz, " cells per group, repl ", i)) + 
                scale_shape_discrete(name = "") + xlab("Cell") + ylab("TPM + 1") + 
                facet_wrap(~plot_group, scales = "free_y") + theme_bw() + 
                theme(legend.position = "bottom", 
                      axis.text.x = element_blank()))
      }, error = function(e) e)
    }
  }
}

char_ds_m <- data.frame(ds = names(char_ds), 
                        n_genes = sapply(char_ds, function(w) w["n_genes"]),
                        silhouette_avg = sapply(char_ds, function(w) w["silhouette_avg"]), 
                        stringsAsFactors = FALSE) %>%
  tidyr::separate(ds, into = c("n_cells", "repl"), sep = "\\.", remove = FALSE) %>%
  dplyr::mutate(n_cells = as.numeric(n_cells)) %>%
  dplyr::mutate(repl = as.numeric(repl)) %>%
  dplyr::arrange(n_cells, repl) %>%
  dplyr::mutate(ds = factor(ds, levels = ds)) %>%
  dplyr::mutate(n_cells = factor(n_cells, levels = unique(n_cells)))
print(char_ds_m %>% ggplot(aes(x = ds, y = n_genes, fill = n_cells)) + geom_bar(stat = "identity") + 
        theme_bw() + xlab("Data set") + ylab("Number of genes") + 
        scale_fill_discrete(name = "Number of cells") + 
        stat_summary(fun.data = function(x) {return(c(y = x, label = x))}, 
                     geom = "text", alpha = 1, size = 2, vjust = -1, 
                     position = position_dodge(width = 0.75)))
print(char_ds_m %>% ggplot(aes(x = ds, y = round(silhouette_avg, 4), fill = n_cells)) + 
        geom_bar(stat = "identity") + 
        theme_bw() + xlab("Data set") + ylab("Average silhouette width") + 
        scale_fill_discrete(name = "Number of cells") + 
        stat_summary(fun.data = function(x) {return(c(y = x, label = x))}, 
                     geom = "text", alpha = 1, size = 2, vjust = -1, 
                     position = position_dodge(width = 0.75)))

char_gene <- Reduce(function(...) dplyr::full_join(..., by = "gene"), char_gene)
char_gene_m <- reshape2::melt(char_gene) %>% 
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group"))) 

char_gene_s <- char_gene_m %>% tidyr::spread(key = mtype, value = value) %>%
  dplyr::mutate(replicate = paste0(ncells, ", repl ", repl)) %>%
  dplyr::arrange(ncells, as.numeric(repl))

char_cells <- Reduce(function(...) dplyr::full_join(..., by = "cell"), char_cells)
char_cells_m <- reshape2::melt(char_cells[, !(colnames(char_cells) %in% 
                                                grep("condition", colnames(char_cells), value = TRUE))]) %>% 
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group")))
char_cells_cond <- reshape2::melt(char_cells[, colnames(char_cells) %in% 
                                               c("cell", grep("condition", 
                                                              colnames(char_cells), value = TRUE))],
                                  id.vars = "cell") %>%
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group"))) %>%
  dplyr::rename(condition = value) %>% dplyr::select(-mtype)
char_cells_m <- dplyr::full_join(char_cells_m, char_cells_cond, by = c("cell", "ncells", "repl"))
char_cells_s <- char_cells_m %>% tidyr::spread(key = mtype, value = value) %>%
  dplyr::mutate(replicate = paste0(ncells, ", repl ", repl)) %>%
  dplyr::arrange(ncells, as.numeric(repl))

                                                
print(ggplot(char_gene_s, aes(x = avecount, y = fraczero)) + geom_point(size = 0.3) +
        theme_bw() + scale_x_log10() + facet_wrap(~forcats::as_factor(replicate)) +
        xlab("Average count per gene") + ylab("Fraction zeros per gene"))
print(ggplot(char_cells_s %>% dplyr::filter(!is.na(condition)), 
             aes(x = libsize, y = fraczero, color = condition)) + geom_point(size = 0.7) +
        theme_bw() + scale_x_log10() + facet_wrap(~forcats::as_factor(replicate)) +
        xlab("Library size") + ylab("Fraction zeros per cell"))

for (tp in c("vartpm", "avecount", "avetpm", "avecensuscount", "avecpm", "varcpm")) {
  nn <- switch(tp,
               vartpm = "Variance of TPM values per gene",
               avecount = "Average count per gene",
               avecensuscount = "Average census count per gene",
               avetpm = "Average TPM per gene",
               avecpm = "Average CPM per gene",
               varcpm = "Variance of CPM values per gene")
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + scale_x_log10() + facet_wrap(~ncells) + 
          theme_bw() + xlab(nn))
}
for (tp in c("fraczero", "fraczerodiff", "fraczerocensus", "cvtpm", "fraczeroround", 
             "cvcpm", "fracimputedup", "fracimputeddown", "fracnonimputed")) {
  nn <- switch(tp, 
               fraczero = "Fraction zeros per gene",
               fraczeroround = "Fraction zeros per gene after rounding",
               fraczerodiff = "Difference (between conditions) of zero fraction per gene",
               cvtpm = "Coefficient of variation (TPM)",
               cvcpm = "Coefficient of variation (CPM)",
               fraczerocensus = "Fraction zeros per gene, census counts",
               fracimputedup = "Fraction of values imputed upwards",
               fracimputeddown = "Fraction of values imputed downwards",
               fracnonimputed = "Fraction of values not imputed")
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + facet_wrap(~ncells) + 
          theme_bw() + xlab(nn))
}
for (tp in c("libsize", "fraczero", "fraczeroround", "libsizecensus", "fraczerocensus")) {
  nn <- switch(tp,
               fraczero = "Fraction zeros per cell",
               fraczeroround = "Fraction zeros per cell after rounding",
               libsize = "Library size per cell",
               libsizecensus = "Library size per cell, census counts",
               fraczerocensus = "Fraction zeros per cell, census counts")
  print(char_cells_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + facet_wrap(~ncells) + 
          theme_bw() + xlab(nn))
  print(char_cells_m %>% dplyr::filter(mtype == tp) %>%
          dplyr::filter(!is.na(value)) %>%
          ggplot(aes(x = condition, y = value)) + geom_boxplot(outlier.size = -1) +
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab(nn) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          facet_wrap(~interaction(ncells, repl)))
  print(char_cells_m %>% dplyr::filter(mtype == tp) %>%
          dplyr::filter(!is.na(value)) %>%
          ggplot(aes(x = cell, y = value, fill = condition)) + geom_bar(stat = "identity") +
          theme_bw() + xlab("Cell") + ylab(nn) + 
          theme(axis.text.x = element_blank(),
                legend.position = "bottom") + 
          facet_wrap(~interaction(ncells, repl)))
}

dev.off()

saveRDS(plots, file = paste0(figdir, "/", dataset, exts, "_dataset_characteristics_plots.rds"))
saveRDS(lapply(list(char_cells_m = char_cells_m, 
                    char_gene_m = char_gene_m,
                    char_ds_m = char_ds_m),
               function(L) {
                 L$dataset <- dataset
                 L$filt <- filt
                 L
               }), 
        paste0(figdir, "/", dataset, exts, "_dataset_characteristics_summary_data.rds"))

sessionInfo()
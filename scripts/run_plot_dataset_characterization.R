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
source("/home/Shared/data/seq/conquer/comparison/scripts/prepare_mae.R")
source("/home/Shared/data/seq/conquer/comparison/scripts/calculate_gene_characteristics.R")
source("/home/Shared/data/seq/conquer/comparison/scripts/calculate_cell_characteristics.R")

if (filt == "") { 
  exts <- filt
} else {
  exts <- paste0("_", filt)
}

pdf(paste0("figures/dataset_characteristics/", dataset, exts, ".pdf"), width = 14, height = 9)

print(dataset)
print(config_file)
print(filt)
print(cell_cycle_file)

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
} else {
  cell_cycle_ids <- cell_cycle_file[[gsub(" ", "_", metadata(mae)$organism)]]
}

sizes <- names(keep_samples)
char_gene <- list()
char_cells <- list()
char_ds <- list()
for (sz in sizes) {
  for (i in 1:nrow(keep_samples[[as.character(sz)]])) {
    message(sz, ".", i)
    L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, filt = filt)
    # if (as.numeric(sz) == max(as.numeric(sizes))) {
    #   ## Save condition information for all cells
    #   cell_group <- L$condt
    # }
    
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
    
    char_ds[[paste0(sz, ".", i)]] <- c(n_genes = nrow(L$count))
    
    ## Plot expression of cell cycle genes
    if (!is.null(cell_cycle_ids)) {
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
    }
  }
}

char_ds_m <- data.frame(ds = names(char_ds), n_genes = sapply(char_ds, function(w) w["n_genes"]),
                        stringsAsFactors = FALSE) %>%
  tidyr::separate(ds, into = c("n_cells", "repl"), sep = "\\.", remove = FALSE) %>%
  dplyr::mutate(n_cells = as.numeric(n_cells)) %>%
  dplyr::mutate(repl = as.numeric(repl)) %>%
  dplyr::arrange(n_cells, repl) %>%
  dplyr::mutate(ds = factor(ds, levels = ds)) %>%
  dplyr::mutate(n_cells = factor(n_cells, levels = unique(n_cells)))
print(char_ds_m %>% ggplot(aes(x = ds, y = n_genes, fill = n_cells)) + geom_bar(stat = "identity") + 
        theme_bw() + xlab("Data set") + ylab("Number of genes") + 
        scale_fill_discrete(name = "Number of cells"))

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
                                               c("cell", grep("condition", colnames(char_cells), value = TRUE))],
                                  id.vars = "cell") %>%
  tidyr::separate(variable, into = c("mtype", "ncells", "repl"), sep = "\\.") %>%
  dplyr::mutate(ncells = factor(paste0(ncells, " cells per group"), 
                                levels = paste0(as.character(sort(as.numeric(unique(ncells)))),
                                                " cells per group"))) %>%
  dplyr::rename(condition = value) %>% dplyr::select(-mtype)
char_cells_m <- dplyr::full_join(char_cells_m, char_cells_cond, by = c("cell", "ncells", "repl"))
                                                
print(ggplot(char_gene_s, aes(x = avecount, y = fraczero)) + geom_point(size = 0.3) +
        theme_bw() + scale_x_log10() + facet_wrap(~forcats::as_factor(replicate)) +
        xlab("Average count per gene") + ylab("Fraction zeros per gene"))

for (tp in c("vartpm", "avecount", "avetpm", "avecensuscount")) {
  nn <- switch(tp,
               vartpm = "Variance of TPM values per gene",
               avecount = "Average count per gene",
               avecensuscount = "Average census count per gene",
               avetpm = "Average TPM per gene")
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + scale_x_log10() + facet_wrap(~ncells) + 
          theme_bw() + xlab(nn))
}
for (tp in c("fraczero", "fraczerodiff", "fraczerocensus", "cvtpm")) {
  nn <- switch(tp, 
               fraczero = "Fraction zeros per gene",
               fraczerodiff = "Difference (between conditions) of zero fraction per gene",
               cvtpm = "Coefficient of variation (TPM)",
               fraczerocensus = "Fraction zeros per gene, census counts")
  print(char_gene_m %>% dplyr::filter(mtype == tp) %>% 
          ggplot(aes(x = value, group = repl, col = repl)) +
          scale_color_discrete(guide = FALSE) + 
          geom_density() + facet_wrap(~ncells) + 
          theme_bw() + xlab(nn))
}
for (tp in c("libsize", "fraczero", "libsizecensus", "fraczerocensus")) {
  nn <- switch(tp,
               fraczero = "Fraction zeros per cell",
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

saveRDS(NULL, paste0("figures/dataset_characteristics/", dataset, exts, ".rds"))

sessionInfo()
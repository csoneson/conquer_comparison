args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(figdir)

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
source("scripts/prepare_mae.R")
source(paste0("scripts/apply_voomlimma.R"))

## First, an example (GSE48968-GPL13112mock, 24 cells per group, repl 1)

config_file <- "config/GSE48968-GPL13112mock.json"
config <- fromJSON(file = config_file)
mae <- readRDS(config$mae)
groupid <- config$groupid
mae <- clean_mae(mae = mae, groupid = groupid)
subsets <- readRDS(config$subfile)
keep_samples <- subsets$keep_samples
imposed_condition <- subsets$out_condition
sz <- "24"
i <- 1

plots <- list()
for (flt in c("", "count_15_25p")) {
  L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, flt)
  dge <- DGEList(L$count, group = L$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~L$condt)
  vm <- voom(dge, design = design, plot = FALSE, save.plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  stopifnot(abs(cor(fit$Amean, vm$voom.xy$x) - 1) <= 1e-15)
  df <- data.frame(Amean = vm$voom.xy$x, 
                   logfc = fit$coefficients[, ncol(fit$coefficients)])
  df2 <- data.frame(logcount = vm$voom.xy$x,
                    sqrtsd = vm$voom.xy$y)
  df3 <- data.frame(logcount = vm$voom.line$x,
                    sqrtsd = vm$voom.line$y)
  
  p1 <- ggplot(df2, aes(x = logcount, y = sqrtsd)) + geom_point(size = 0.5) + 
    geom_line(data = df3, col = "red", size = 2) + 
    xlab("log2(count + 0.5)") + ylab("sqrt(standard deviation)") + theme_bw() + 
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  p2 <- ggplot(df, aes(x = Amean, y = logfc)) + 
    geom_point(size = 0.5) + geom_hline(yintercept = 0, color = "black") + 
    geom_smooth(se = FALSE, size = 2) + xlab("log2(count + 0.5)") + 
    ylab("log-fold-change") + theme_bw() + 
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  if (flt == "") labs <- c("A", "B")
  else labs <- c("C", "D")
  plots[[paste0("flt_", flt)]] <- plot_grid(p1, p2, rel_widths = c(1, 1), labels = labs)
}

titles <- list(
  flt_ = ggdraw() + 
    draw_label("GSE48968-GPL13112null, 24 cells per group, no filtering",
               fontface = "bold"),
  flt_count_15_25p = ggdraw() + 
    draw_label("GSE48968-GPL13112null, 24 cells per group, count>15 in >25% of cells",
               fontface = "bold")
)

## -------------------------------------------------------------------------- ##
## Next, go through all data sets. Record the FPR and the fraction of genes to 
## the left of the peak of the smoothing line in the voom plot.
datasets <- structure(c("GSE45719mock", "GSE74596mock", "EMTAB2805mock", 
                        "GSE60749-GPL13112mock", "GSE48968-GPL13112mock",
                        "UsoskinGSE59739mock"), 
                      names = c("GSE45719mock", "GSE74596mock", "EMTAB2805mock", 
                                "GSE60749-GPL13112mock", "GSE48968-GPL13112mock", 
                                "UsoskinGSE59739mock"))
FPR <- c()
frac_below_peak <- c()
peak_height <- c()
peak_dist_to_min <- c()
peak_rel_dist_to_min <- c()
for (ds in datasets) {
  message(ds)
  config_file <- paste0("config/", ds, ".json")
  config <- fromJSON(file = config_file)
  mae <- readRDS(config$mae)
  groupid <- config$groupid
  mae <- clean_mae(mae = mae, groupid = groupid)
  subsets <- readRDS(config$subfile)
  keep_samples <- subsets$keep_samples
  imposed_condition <- subsets$out_condition
  for (flt in c("", "TPM_1_25p", "TPM_5_25p", "TPM_11_25p")) {
    for (sz in names(imposed_condition)) {
      for (i in 1:nrow(imposed_condition[[sz]])) {
        L <- subset_mae(mae, keep_samples, sz, i, imposed_condition, flt)
        dge <- DGEList(L$count, group = L$condt)
        dge <- calcNormFactors(dge)
        design <- model.matrix(~L$condt)
        vm <- voom(dge, design = design, plot = FALSE, save.plot = TRUE)
        fit <- lmFit(vm, design = design)
        fit <- eBayes(fit)
        tt <- topTable(fit, n = Inf, adjust.method = "BH")
        FPR[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- mean(tt$P.Value <= 0.05)
        frac_below_peak[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- 
          mean(vm$voom.line$x <= vm$voom.line$x[which.max(vm$voom.line$y)])
        peak_height[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- 
          max(vm$voom.line$y) - vm$voom.line$y[1]
        peak_dist_to_min[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- 
          vm$voom.line$x[which.max(vm$voom.line$y)] - min(vm$voom.line$x)
        peak_rel_dist_to_min[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- 
          (vm$voom.line$x[which.max(vm$voom.line$y)] - min(vm$voom.line$x))/
          (max(vm$voom.line$x) - min(vm$voom.line$x))
      }
    }
  }  
}
df4 <- Reduce(function(...) dplyr::full_join(..., by = "nm"), 
              list(data.frame(nm = names(peak_height), peak_height),
                   data.frame(nm = names(FPR), FPR),
                   data.frame(nm = names(frac_below_peak), frac_below_peak),
                   data.frame(nm = names(peak_dist_to_min), peak_dist_to_min),
                   data.frame(nm = names(peak_rel_dist_to_min), peak_rel_dist_to_min))) %>%
  tidyr::separate(nm, into = c("dataset", "filt", "ncells", "repl"), sep = "_") %>%
  dplyr::mutate(dataset = gsub("mock$", "null", dataset))

p41 <- ggplot(df4, aes(x = peak_height, y = FPR, shape = dataset)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + geom_point() + 
  geom_smooth() + theme_bw() + xlab("Peak height in the voom plot")
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))
p42 <- ggplot(df4, aes(x = frac_below_peak, y = FPR, shape = dataset)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + geom_point() + 
  geom_smooth() + theme_bw() + 
  xlab("Fraction of genes to the left of the peak in the voom plot") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))
p43 <- ggplot(df4, aes(x = frac_below_peak, y = FPR, shape = dataset)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + geom_point() + 
  geom_smooth() + ylim(0, 0.2) + theme_bw() + 
  xlab("Fraction of genes to the left of the peak in the voom plot") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))
p44 <- ggplot(df4, aes(x = peak_dist_to_min, y = FPR, shape = dataset)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + geom_point() + 
  geom_smooth() + theme_bw() + 
  xlab("Distance between peak and leftmost point in the voom plot") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))
p45 <- ggplot(df4, aes(x = peak_rel_dist_to_min, y = FPR, shape = dataset)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + geom_point() + 
  geom_smooth() + theme_bw() + 
  xlab("Relative position of the peak between the leftmost and rightmost points in the voom plot") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

pdf(paste0(figdir, "/voomlimma_investigation.pdf"), width = 12, height = 18)
print(plot_grid(plot_grid(titles[["flt_"]], plots[["flt_"]], ncol = 1, rel_heights = c(0.1, 1)),
                plot_grid(titles[["flt_count_15_25p"]], plots[["flt_count_15_25p"]], 
                          ncol = 1, rel_heights = c(0.1, 1)),
                p42, ncol = 1, rel_heights = c(1, 1), labels = c("", "", "E")))
dev.off()

saveRDS(NULL, file = paste0(figdir, "/voomlimma_investigation.rds"))
sessionInfo()
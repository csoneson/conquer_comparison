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
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
source("scripts/prepare_mae.R")

## First, an example (GSE48968-GPL13112mock, 24 cells per group, repl 1)

config <- fromJSON(file = "config/GSE48968-GPL13112mock.json")
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
  ## Subset data set and apply voom/limma
  L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i,
                  imposed_condition = imposed_condition, filt = flt, 
                  impute = config$impute)
  dge <- DGEList(L$count, group = L$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~L$condt)
  vm <- voom(dge, design = design, plot = FALSE, save.plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  ## Prepare plots
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
    geom_smooth(se = FALSE, size = 2, method = "loess") + xlab("log2(count + 0.5)") + 
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
                        "UsoskinGSE59739mock", "GSE62270-GPL17021mock", 
                        "10XMonoCytoTmock"), 
                      names = c("GSE45719mock", "GSE74596mock", "EMTAB2805mock", 
                                "GSE60749-GPL13112mock", "GSE48968-GPL13112mock", 
                                "UsoskinGSE59739mock", "GSE62270-GPL17021mock", 
                                "10XMonoCytoTmock"))
FPR <- c()
frac_below_peak <- c()
for (ds in datasets) {
  message(ds)
  config <- fromJSON(file = paste0("config/", ds, ".json"))
  mae <- readRDS(config$mae)
  groupid <- config$groupid
  mae <- clean_mae(mae = mae, groupid = groupid)
  subsets <- readRDS(config$subfile)
  keep_samples <- subsets$keep_samples
  imposed_condition <- subsets$out_condition
  for (flt in c("", "TPM_1_25p", "TPM_5_25p", "TPM_11_25p")) {
    for (sz in names(imposed_condition)) {
      for (i in 1:nrow(imposed_condition[[sz]])) {
        L <- subset_mae(mae = mae, keep_samples = keep_samples, sz = sz, i = i, 
                        imposed_condition = imposed_condition, filt = flt, 
                        impute = config$impute)
        dge <- DGEList(L$count, group = L$condt)
        dge <- calcNormFactors(dge)
        design <- model.matrix(~L$condt)
        vm <- voom(dge, design = design, plot = FALSE, save.plot = TRUE)
        fit <- lmFit(vm, design = design)
        fit <- eBayes(fit)
        tt <- topTable(fit, n = Inf, adjust.method = "BH")
        FPR[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- mean(tt$P.Value < 0.05)
        frac_below_peak[paste0(ds, "_", gsub("_", ".", flt), "_", sz, "_", i)] <- 
          mean(vm$voom.line$x <= vm$voom.line$x[which.max(vm$voom.line$y)])
      }
    }
  }  
}
df4 <- Reduce(function(...) dplyr::full_join(..., by = "nm"), 
              list(data.frame(nm = names(FPR), FPR),
                   data.frame(nm = names(frac_below_peak), frac_below_peak))) %>%
  tidyr::separate(nm, into = c("dataset", "filt", "ncells", "repl"), sep = "_") %>%
  dplyr::mutate(dataset = gsub("mock$", "null", dataset))

pch <- c(16, 17, 15, 3, 7, 8, 4, 6, 9, 10, 11, 12, 13, 14, 1, 2, 5, 18, 19, 20)[1:length(unique(df4$dataset))]
names(pch) <- unique(df4$dataset)

p3 <- ggplot(df4, aes(x = frac_below_peak, y = FPR)) + 
  geom_hline(yintercept = 0.05, color = "red", size = 2) + 
  geom_point(aes(shape = dataset)) + 
  scale_shape_manual(values = pch, name = "Data set") + 
  geom_smooth(method = "loess") + theme_bw() + 
  xlab("Fraction of genes to the left of the peak in the voom plot") + 
  ylab("FPR (fraction of genes with p < 0.05)") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

pdf(paste0(figdir, "/voomlimma_investigation.pdf"), width = 12, height = 18)
print(plot_grid(plot_grid(titles[["flt_"]], plots[["flt_"]], ncol = 1, rel_heights = c(0.1, 1)),
                plot_grid(titles[["flt_count_15_25p"]], plots[["flt_count_15_25p"]], 
                          ncol = 1, rel_heights = c(0.1, 1)),
                p3, ncol = 1, rel_heights = c(1, 1, 1), labels = c("", "", "E")))
dev.off()

saveRDS(NULL, file = paste0(figdir, "/voomlimma_investigation.rds"))
sessionInfo()
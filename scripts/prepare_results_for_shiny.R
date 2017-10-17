args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(fracnards)
print(nbrgenesrds)
print(type1errorrds)
print(trueperfrds)
print(timingrds)
print(decharacrds)
print(origvsmockrds)
print(perfsummaryrds)
print(crossmethodconsrds)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(forcats))
source("scripts/plot_setup.R")

ttest <- function(x, y) {
  tryCatch({
    t.test(x, y, var.equal = FALSE)$stat
  }, error = function(e) NA)
}

fracnards <- readRDS(fracnards) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(fracNA = round(fracNA, 3)) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, fracNA, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact)
nbrgenesrds <- readRDS(nbrgenesrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(nbr_sign_adjp0.05_rel = round(nbr_sign_adjp0.05_rel, 3)) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, nbr_sign_adjp0.05_rel, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact)
type1errorrds <- readRDS(type1errorrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(FPR = round(FPR, 3)) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, FPR, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact) %>%
  dplyr::mutate(dataset = gsub("mock", "null", dataset))
fdprds <- readRDS(trueperfrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::rename(FDP = FDR) %>%
  dplyr::select(-AUROC) %>%
  dplyr::mutate(FDP = round(FDP, 3)) %>%
  dplyr::mutate(TPR = round(TPR, 3)) %>%
  dplyr::mutate(fracNA = round(fracNA, 3)) %>% 
  dplyr::mutate(method = forcats::fct_reorder(method, FDP, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact)
tprrds <- readRDS(trueperfrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::rename(FDP = FDR) %>%
  dplyr::select(-AUROC) %>%
  dplyr::mutate(TPR = round(TPR, 3)) %>%
  dplyr::mutate(FDP = round(FDP, 3)) %>%
  dplyr::mutate(fracNA = round(fracNA, 3)) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, TPR, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact)
aurocrds <- readRDS(trueperfrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(AUROC = round(AUROC, 3)) %>%
  dplyr::mutate(fracNA = round(fracNA, 3)) %>%
  dplyr::select(-FDR, -TPR) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, AUROC, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact)
timingrds <- readRDS(timingrds)$timing %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(rel_timing = round(rel_timing, 5)) %>%
  dplyr::select(-ngenes_cat) %>%
  dplyr::mutate(method = forcats::fct_reorder(method, rel_timing, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(ncells = ncells_fact) %>%
  dplyr::mutate(dataset = gsub("mock", "null", dataset))
decharacrds <- readRDS(decharacrds) %>%
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(snr = round(snr, 3)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::rename(characteristic = charac)
origvsmockrds <- readRDS(origvsmockrds) %>% 
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::mutate(AUCs = round(AUCs, 3)) %>%
  dplyr::filter(k == 100) %>%
  dplyr::group_by(method, dataset, filt, ncells, k) %>%
  dplyr::mutate(tokeep = length(unique(tp)) == 2) %>%
  dplyr::filter(tokeep) %>% 
  dplyr::mutate(signal_vs_mock = mean(AUCs[tp == "signal"]) - mean(AUCs[tp == "mock"])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = forcats::fct_reorder(method, signal_vs_mock, 
                                              fun = median, na.rm = TRUE,
                                              .desc = TRUE)) %>%
  dplyr::rename(filtering = filt) %>%
  dplyr::mutate(tp = replace(tp, tp == "mock", "null")) %>%
  dplyr::filter(filtering == "TPM_1_25p") %>%
  dplyr::mutate(repl = paste0(replicate1, ".", replicate2)) %>%
  dplyr::select(-replicate1, -replicate2, -tokeep) %>%
  dplyr::rename(AUCC = AUCs) %>%
  dplyr::mutate(ncells = droplevels(ncells_fact)) %>%
  dplyr::rename(signaltype = tp) %>%
  dplyr::mutate(signaltype = factor(signaltype, levels = c("signal", "null")))
crossmethodconsrds <- readRDS(crossmethodconsrds) %>% 
  dplyr::filter(k == 100) %>%
  dplyr::select(dataset, filt, method1, method2, ncells, repl, k, AUCs) %>%
  dplyr::mutate(filt = replace(filt, filt == "", "No filtering")) %>%
  dplyr::rename(filtering = filt) %>% dplyr::ungroup()
cm1 <- crossmethodconsrds
cm2 <- crossmethodconsrds
cm2[, c("method1", "method2")] <- cm2[, c("method2", "method1")]
crossmethodconsrds <- rbind(cm1, cm2)
  
perfsummaryrds <- readRDS(perfsummaryrds)

allmethods <- sort(unique(c(as.character(fracnards$method), as.character(nbrgenesrds$method), 
                            as.character(type1errorrds$method), 
                            as.character(fdprds$method), as.character(tprrds$method), 
                            as.character(aurocrds$method), as.character(timingrds$method),
                            as.character(decharacrds$method), as.character(origvsmockrds$method),
                            as.character(crossmethodconsrds$method1))))
allncells <- as.character(sort(as.numeric(
  unique(c(as.character(fracnards$ncells), as.character(nbrgenesrds$ncells), 
           as.character(type1errorrds$ncells), 
           as.character(fdprds$ncells), as.character(tprrds$ncells), 
           as.character(aurocrds$ncells), as.character(timingrds$ncells),
           as.character(decharacrds$ncells), as.character(origvsmockrds$ncells),
           as.character(crossmethodconsrds$ncells))))))
alldatasets <- sort(unique(c(fracnards$dataset, nbrgenesrds$dataset, type1errorrds$dataset, 
                             fdprds$dataset, tprrds$dataset, aurocrds$dataset, 
                             timingrds$dataset, decharacrds$dataset, origvsmockrds$dataset,
                             crossmethodconsrds$dataset)))
allfilterings <- sort(unique(c(fracnards$filtering, nbrgenesrds$filtering, type1errorrds$filtering, 
                               fdprds$filtering, tprrds$filtering, aurocrds$filtering, 
                               timingrds$filtering, decharacrds$filtering, origvsmockrds$filtering,
                               crossmethodconsrds$filtering)))

saveRDS(list(fracna = list(id = "Fraction NA adj.p", 
                           data = fracnards,
                           yvar = "fracNA", 
                           ymin = 0, ymax = 1,
                           facet2 = c("filtering", "none")),
             nbrgenes = list(id = "Number of DEGs", 
                             data = nbrgenesrds,
                             yvar = "nbr_sign_adjp0.05",
                             ymin = 0, ymax = max(nbrgenesrds$nbr_sign_adjp0.05),
                             facet2 = c("filtering", "none")),
             type1error = list(id = "Type 1 error",
                               data = type1errorrds,
                               yvar = "FPR", 
                               ymin = 0, ymax = 1,
                               facet2 = c("filtering", "none")),
             fdp = list(id = "FDP at adj.p = 0.05 cutoff",
                        data = fdprds,
                        yvar = "FDP",
                        ymin = 0, ymax = 1,
                        facet2 = c("filtering", "none")),
             tpr = list(id = "TPR at adj.p = 0.05 cutoff",
                        data = tprrds,
                        yvar = "TPR",
                        ymin = 0, ymax = 1,
                        facet2 = c("filtering", "none")),
             auroc = list(id = "AUROC",
                          data = aurocrds,
                          yvar = "AUROC",
                          ymin = 0, ymax = 1,
                          facet2 = c("filtering", "none")),
             timing = list(id = "Time requirement",
                           data = timingrds,
                           yvar = "rel_timing",
                           ymin = 0, ymax = 1,
                           facet2 = c("filtering", "none")),
             decharac = list(id = "DE genes characteristics",
                             data = decharacrds,
                             yvar = "snr",
                             ymin = min(decharacrds$snr), ymax = max(decharacrds$snr)),
             origvsmock = list(id = "Concordance scores",
                               data = origvsmockrds,
                               yvar = "AUCC",
                               ymin = 0, ymax = 1,
                               facet2 = c("signaltype")), 
             crossmethodcons = crossmethodconsrds, 
             perfsummary = perfsummaryrds, 
             cols = cols, allmethods = allmethods, allfilterings = allfilterings, 
             allncells = allncells, alldatasets = alldatasets), file = outrds)

sessionInfo()
date()

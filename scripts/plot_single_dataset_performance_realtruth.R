source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

plot_performance_realtruth <- function(cobra, truth, colvec, exts = exts, summary_data = list()) {
  cobra <- COBRAData(truth = truth, 
                     object_to_extend = cobra)
  cobraperf <- calculate_performance(cobra, binary_truth = "status", 
                                     aspects = c("fdrtpr", "fdrtprcurve", "fpr", "tpr"), thrs = c(0.05))
  
  sizes_ncells <- unique(paste(get_nsamples(basemethods(cobraperf)), 
                               get_repl(basemethods(cobraperf)), sep = "."))
  for (szi in sizes_ncells) {
    spl <- strsplit(szi, "\\.")[[1]]
    keepmethods <- basemethods(cobraperf)[get_nsamples(basemethods(cobraperf)) == spl[1] 
                                          & get_repl(basemethods(cobraperf)) == spl[2]]
    c2 <- structure(colvec, names = paste0(names(colvec),
                                           ".", spl[1], ".", spl[2]))
    cobraplot <- 
      prepare_data_for_plot(cobraperf, 
                            colorscheme = c2[keepmethods],
                            keepmethods = keepmethods)
    print(plot_fdrtprcurve(cobraplot, plottype = c("point", "curve")))
  }

  print(fdrtpr(cobraperf) %>% tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.") %>%
          ggplot(aes(x = method, y = FDR, color = method, shape = ncells)) + 
          geom_hline(yintercept = 0.05) + 
          geom_point() + theme_bw() + 
          scale_color_manual(values = structure(colvec, names = gsub(exts, "", names(colvec))), name = "") + 
          theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 13)) + 
          xlab(""))
  return(invisible(summary_data))
}

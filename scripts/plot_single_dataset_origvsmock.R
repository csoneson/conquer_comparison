plot_compare_orig_mock <- function(concordances, colvec, K0, summary_data = list()) {
  concs <- lapply(names(concordances), function(ncbr) {
    concordances[[ncbr]]$concordance_pairwise_bymethod %>%
      dplyr::filter(k %in% K0) %>%
      dplyr::mutate(tp = ncbr) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_mock", "mock")) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_", "signal"))
  })
  concs <- do.call(rbind, concs)
  
  summary_data$concordances <- concs
  
  for (k0 in K0) {
    for (nbrsamples in unique(intersect(subset(concs, tp == "signal")$ncells,
                                        subset(concs, tp == "mock")$ncells))) {
      message(nbrsamples)
      print(concs %>% dplyr::filter(k == k0) %>% dplyr::filter(ncells == nbrsamples) %>%
              ggplot(aes(x = method, y = AUCs, color = method, shape = tp)) + 
              geom_point(size = 5) + theme_bw() + xlab("") + 
              ylab(paste0("Area under concordance curve between data set instances, top-", k0, " genes")) + 
              scale_color_manual(values = colvec, name = "") + 
              scale_shape_discrete(name = "") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              ggtitle(paste0(nbrsamples, " cells per group")))
    }
  }
  
  return(invisible(summary_data))
}

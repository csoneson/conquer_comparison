source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

plot_compare_orig_mock <- function(concordances, colvec, summary_data = list()) {
  concs <- lapply(names(concordances), function(ncbr) {
    conc <- concordances[[ncbr]]$concordance_pairwise_auc
    conc %>%
      dplyr::mutate(tp = ncbr) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_mock", "mock")) %>%
      dplyr::mutate(tp = replace(tp, tp == "tp_", "original"))
  })
  concs <- do.call(rbind, concs)
  
  summary_data$concordances <- concs
  
  for (nbrsamples in unique(intersect(subset(concs, tp == "original")$ncells1,
                                      subset(concs, tp == "mock")$ncells1))) {
    print(concs %>% dplyr::filter(ncells1 == nbrsamples & ncells2 == nbrsamples) %>%
            ggplot(aes(x = method, y = auc, color = method, shape = tp)) + 
            geom_point(size = 5) + theme_bw() + xlab("") + 
            ylab(paste0("Area under concordance curve")) + 
            scale_color_manual(values = colvec, name = "") + 
            scale_shape_discrete(name = "") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
            ggtitle(paste0(nbrsamples, " cells per group")))
  }
  
  return(invisible(summary_data))
}

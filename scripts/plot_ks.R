source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Plot p-value and statistic for Kolmogorov-Smirnov test for uniformity of
#' p-values
#' 
plot_ks <- function(cobra, colvec, exts = exts, summary_data = list()) {
  ## For each method with p-values, calculate p-value from K-S test for uniformity
  pvs <- pval(cobra)
  ksp <- apply(pvs, 2, function(x) ks.test(x = x[!is.na(x)], y = punif, min = 0, max = 1)$p.value)
  kst <- apply(pvs, 2, function(x) ks.test(x = x[!is.na(x)], y = punif, min = 0, max = 1)$statistic)
  
  print(data.frame(method = names(kst), KSstat = kst, stringsAsFactors = FALSE) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = 
                          factor(nbr_samples, 
                                 levels = as.character(sort(as.numeric(as.character(nbr_samples)))))) %>% 
          ggplot(aes(x = method, y = KSstat, color = method, shape = nbr_samples)) + 
          geom_point(size = 5) + 
          theme_bw() + 
          xlab("") + ylab("Kolmogorov-Smirnov statistic, p-value distribution") + 
          scale_color_manual(values = colvec, name = "") + 
          scale_shape_discrete(name = "Number of cells") + 
          theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 13),
                axis.title.y = element_text(size = 15)))
  
  for (i in sort(unique(as.numeric(get_nsamples(names(kst)))))) {
    print(data.frame(method = names(kst), KSstat = kst, stringsAsFactors = FALSE) %>% 
            tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
            dplyr::filter(nbr_samples == i) %>% 
            ggplot(aes(x = method, y = KSstat, color = method)) + 
            geom_point(size = 5) + 
            theme_bw() + 
            xlab("") + ylab("Kolmogorov-Smirnov statistic, p-value distribution") + 
            scale_color_manual(values = colvec, name = "") + 
            theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 13),
                  axis.title.y = element_text(size = 15)) +
            ggtitle(paste0(i, " cells per condition")))
  }
  return(invisible(summary_data))
}

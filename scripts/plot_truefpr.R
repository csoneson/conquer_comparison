source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Calculate and plot "true FPR" (fraction of genes with nominal p<0.05 in null
#' situation)
#' 
plot_truefpr <- function(cobra, colvec, exts = exts, summary_data = list()) {
  pvs <- pval(cobra)
  fpr <- apply(pvs, 2, function(x) length(which(x < 0.05))/length(x[!is.na(x)]))
  
  summary_data$fracpbelow0.05 <- 
    rbind(summary_data$fracpbelow0.05, 
          data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE))
  
  print(data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE) %>% 
          tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
          dplyr::mutate(nbr_samples = 
                          factor(nbr_samples, 
                                 levels = as.character(sort(as.numeric(as.character(nbr_samples)))))) %>% 
          ggplot(aes(x = method, y = FPR, color = method, shape = nbr_samples)) + geom_point(size = 5) + 
          theme_bw() + xlab("") + ylab("Fraction of genes with nominal p < 0.05") + 
          geom_hline(yintercept = 0.05) + 
          scale_color_manual(values = colvec, name = "") +
          scale_shape_discrete(name = "Number of cells") + 
          theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 13),
                axis.title.y = element_text(size = 15)))
  
  for (i in sort(unique(as.numeric(get_nsamples(names(fpr)))))) {
    print(data.frame(method = names(fpr), FPR = fpr, stringsAsFactors = FALSE) %>% 
            tidyr::separate(method, into = c("method", "nbr_samples", "replicate"), sep = "\\.") %>%
            dplyr::filter(nbr_samples == i) %>% 
            ggplot(aes(x = method, y = FPR, color = method)) + geom_point(size = 5) + 
            theme_bw() + xlab("") + ylab("Fraction of genes with nominal p < 0.05") + 
            geom_hline(yintercept = 0.05) + 
            scale_color_manual(values = colvec, name = "method") + 
            theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 13),
                  axis.title.y = element_text(size = 15)) +
            ggtitle(paste0(i, " cells per condition")))
  }
  return(invisible(summary_data))
}

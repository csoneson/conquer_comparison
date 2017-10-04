plot_truefpr <- function(cobra, colvec, exts, summary_data = list()) {
  pvs <- pval(cobra)
  fpr <- apply(pvs, 2, function(x) length(which(x <= 0.05))/length(x[!is.na(x)]))
  nsign <- apply(pvs, 2, function(x) length(which(x <= 0.05)))
  ntest <- apply(pvs, 2, function(x) length(x[!is.na(x)]))
  
  df <- data.frame(method = names(fpr), FPR = fpr, nsign = nsign, ntest = ntest, 
                   stringsAsFactors = FALSE)
  
  summary_data$fracpbelow0.05 <- df
  
  print(df %>% 
          tidyr::separate(method, into = c("method", "ncells", "replicate"), sep = "\\.") %>%
          dplyr::mutate(method = gsub(exts, "", method)) %>%
          dplyr::mutate(
            ncells = factor(ncells, 
                            levels = as.character(sort(unique(as.numeric(as.character(ncells))))))) %>% 
          ggplot(aes(x = method, y = FPR, color = method, shape = ncells)) + geom_hline(yintercept = 0.05) + 
          geom_point(size = 5) + 
          theme_bw() + xlab("") + ylab("Fraction of tested genes with nominal p < 0.05") + 
          scale_color_manual(values = structure(colvec, names = gsub(exts, "", names(colvec))), name = "") +
          scale_shape_discrete(name = "Number of cells\nper group") + 
          theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 13),
                axis.title.y = element_text(size = 15)))
  
  for (i in sort(unique(as.numeric(get_nsamples(names(fpr)))))) {
    print(df %>% 
            tidyr::separate(method, into = c("method", "ncells", "replicate"), sep = "\\.") %>%
            dplyr::mutate(method = gsub(exts, "", method)) %>%
            dplyr::filter(ncells == i) %>% 
            ggplot(aes(x = method, y = FPR, color = method)) + geom_hline(yintercept = 0.05) + 
            geom_point(size = 5) + 
            theme_bw() + xlab("") + ylab("Fraction of tested genes with nominal p < 0.05") + 
            scale_color_manual(values = structure(colvec, names = gsub(exts, "", names(colvec))), name = "") + 
            theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 13),
                  axis.title.y = element_text(size = 15)) +
            ggtitle(paste0(i, " cells per group")))
  }
  return(invisible(summary_data))
}

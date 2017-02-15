source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

#' Plot distribution of gene characteristics for significant/non-significant
#' genes found with each method
#' 
plot_results_characterization <- function(cobra, colvec, summary_data = list()) {
  sizes_nsamples <- gsub("avetpm.", "", grep("avetpm.", colnames(truth(cobra)), value = TRUE))
  for (szi in sizes_nsamples) {
    message(szi)
    pvals <- pval(cobra)[, grep(paste0("\\.", szi, "$"), colnames(pval(cobra)))]
    pvals[is.na(pvals)] <- 1
    padjs <- padj(cobra)[, grep(paste0("\\.", szi, "$"), colnames(padj(cobra)))]
    padjs[is.na(padjs)] <- 1
    tth <- truth(cobra)[, grep(paste0("\\.", szi, "$"), colnames(truth(cobra)))]
    colnames(tth) <- gsub(paste0("\\.", szi), "", colnames(tth))
    
    ## Number of methods calling each variable significant
    tmp1 <- rowSums(padjs[match(rownames(tth), rownames(padjs)), ] <= 0.05)
    tth$nbr_methods <- tmp1[match(rownames(tth), names(tmp1))]
    
    df2 <- merge(melt(as.matrix(padjs)), tth, by.x = "Var1", by.y = 0, all = TRUE)
    df2 <- subset(df2, rowSums(is.na(df2)) == 0)
    df2$sign <- df2$value <= 0.05
    df2$Var2 <- factor(df2$Var2, levels = sort(levels(df2$Var2)))
    col_reorder <- colvec
    names(col_reorder) <- paste0(names(col_reorder), ".", szi)
    for (y in c("avetpm", "avecount", "vartpm")) {
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~log2(x), x = as.name(y))) %>%
                group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch((median(newy[sign]) - median(newy[!sign]))/median(newy),
                                        error = function(e) NA)) %>%
                dplyr::mutate(charac = paste0("log2_", y)))
      
      print(ggplot(df2, aes_string(x = "Var2", y = paste0("log2(", y, ")"), 
                                   fill = "Var2", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = log2(max(df2[, y])),
                                                            label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              ylab(paste0(expression(log[2]), "(", ifelse(y == "avetpm", "average TPM)", 
                                                          ifelse(y == "avecount", "average count)", 
                                                                 "variance(TPM))")))))
      print(ggplot(df2, aes_string(x = paste0("log2(", y, ")"), fill = "Var2", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var2, scales = "free_y") + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(paste0(expression(log[2]), "(", ifelse(y == "avetpm", "average TPM)", 
                                                          ifelse(y == "avecount", "average count)", 
                                                                 "variance(TPM))")))))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = paste0("log2(", y, ")"), 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + ylab(paste0(expression(log[2]), "(",
                                                        ifelse(y == "avetpm", "average TPM)", 
                                                               ifelse(y == "avecount", "average count)",
                                                                      "variance(TPM))")))) + 
              xlab("Number of methods calling gene significant"))
    }
    for (y in c("fraczero", "fraczerodiff")) {
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~x, x = as.name(y))) %>%
                group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch((median(newy[sign]) - median(newy[!sign]))/median(newy),
                                        error = function(e) NA)) %>%
                dplyr::mutate(charac = y))
      
      print(ggplot(df2, aes_string(x = "Var2", y = y, fill = "Var2", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = max(df2[, y]), label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              ylab(ifelse(y == "fraczero", "Zero fraction", "Difference in zero fraction")))
      print(ggplot(df2, aes_string(x = y, fill = "Var2", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = col_reorder, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var2, scales = "free_y") + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(ifelse(y == "fraczero", "Zero fraction", "Difference in zero fraction")))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = y, 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + ylab(ifelse(y == "fraczero", "Zero fraction", 
                                                        "Difference in zero fraction")) + 
              xlab("Number of methods calling gene significant"))
    }
  }
  return(invisible(summary_data))
}

source("/home/Shared/data/seq/conquer/comparison/scripts/plot_setup.R")

plot_results_characterization <- function(cobra, colvec, exts, summary_data = list()) {
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
    
    df2 <- dplyr::full_join(melt(as.matrix(padjs)), cbind(Var1 = rownames(tth), tth), 
                            by = "Var1")
    df2 <- subset(df2, rowSums(is.na(df2 %>% dplyr::select(-avecensuscount, -fraczerocensus))) == 0)
    df2$sign <- df2$value <= 0.05
    
    df2$Var3 <- gsub(exts, "", gsub(paste0(".", szi), "", df2$Var2))
    names(colvec) <- gsub(exts, "", names(colvec))
    for (y in c("avetpm", "avecount", "vartpm")) {
      nn <- switch(y,
                   avetpm = "average TPM",
                   avecount = "average count",
                   vartpm = "variance(TPM)")
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~log2(x), x = as.name(y))) %>%
                dplyr::group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch(median(newy[sign]) - median(newy[!sign]), error = function(e) NA)) %>%
                dplyr::mutate(charac = paste0("log2_", y)))
      
      print(ggplot(df2, aes_string(x = "Var3", y = paste0("log2(", y, ")"), 
                                   fill = "Var3", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                    axis.text.y = element_text(size = 12),
                    axis.title.y = element_text(size = 13)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = log2(max(df2[, y])),
                                                            label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              xlab("") + 
              ylab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn))) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
      print(ggplot(df2, aes_string(x = paste0("log2(", y, ")"), fill = "Var3", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var3, scales = "free_y") + 
              theme(axis.title.y = element_text(size = 13)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn))) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = paste0("log2(", y, ")"), 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + 
              ylab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn))) + 
              xlab("Number of methods calling gene significant") + 
              theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 13)) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
    }
    for (y in c("fraczero", "fraczeroround", "fraczerodiff", "cvtpm")) {
      nn <- switch(y, 
                   fraczero = "Fraction zeros",
                   fraczeroround = "Fraction zeros after rounding",
                   fraczerodiff = "Difference in fraction zeros between conditions",
                   cvtpm = "coefficient of variation (TPM)")
      ## Populate summary_data with t-statistic between significant and non-significant genes
      summary_data$stats_charac <- 
        rbind(summary_data$stats_charac, 
              df2 %>% dplyr::mutate_(newy = interp(~x, x = as.name(y))) %>%
                group_by(Var2) %>% 
                dplyr::summarise(
                  tstat = tryCatch(t.test(newy[sign], newy[!sign])$statistic, error = function(e) NA),
                  mediandiff = tryCatch(median(newy[sign]) - median(newy[!sign]), error = function(e) NA)) %>%
                dplyr::mutate(charac = y))
      
      print(ggplot(df2, aes_string(x = "Var3", y = y, fill = "Var3", dodge = "sign", alpha = "sign")) + 
              geom_boxplot() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                    axis.text.y = element_text(size = 12),
                    axis.title.y = element_text(size = 13)) + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              stat_summary(fun.data = function(x) {return(c(y = max(df2[, y]), label = length(x)))}, 
                           geom = "text", alpha = 1, size = 2, vjust = -1, 
                           position = position_dodge(width = 0.75)) + 
              xlab("") + ylab(nn) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
      print(ggplot(df2, aes_string(x = y, fill = "Var3", alpha = "sign")) + 
              geom_density() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
              scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
              facet_wrap(~Var3, scales = "free_y") + 
              guides(alpha = guide_legend(override.aes = 
                                            list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                                 colour = NA))) + 
              xlab(nn) + 
              theme(axis.title.y = element_text(size = 13)) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
      print(ggplot(tth, aes_string(x = "nbr_methods", y = y, 
                                   group = "nbr_methods")) + 
              geom_boxplot() + theme_bw() + ylab(nn) + 
              xlab("Number of methods calling gene significant") + 
              theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 13)) + 
              ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi)))
    }
  }
  return(invisible(summary_data))
}

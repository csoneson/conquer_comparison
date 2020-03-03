tstatmod <- function(x, y, lmin) {
  ## Perform t-test only if x and y both have more than lmin non-NA elements
  if (length(x[!is.na(x)]) >= lmin & length(y[!is.na(y)]) >= lmin) t.test(x, y)$statistic
  else NA
}

snr <- function(x, y, lmin) {
  ## Calculate SNR if x and y both have more than lmin non-NA elements
  if (length(x[!is.na(x)]) >= lmin & length(y[!is.na(y)]) >= lmin) {
    (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))/(sd(x, na.rm = TRUE) + sd(y, na.rm = TRUE))
  } else {
    NA
  }
}

mediandiff <- function(x, y, lmin) {
  ## Calculate the difference between medians if x and y both have more than lmin non-NA elements
  if (length(x[!is.na(x)]) >= lmin & length(y[!is.na(y)]) >= lmin) 
    median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
  else NA
}

plot_results_characterization <- function(cobra, colvec, exts, summary_data = list()) {
  sizes_nsamples <- gsub("avetpm.", "", grep("avetpm.", colnames(truth(cobra)), value = TRUE))
  
  for (szi in sizes_nsamples) {
    message(szi)
    pvals <- pval(cobra)[, grep(paste0("\\.", szi, "$"), colnames(pval(cobra)))]
    padjs <- padj(cobra)[, grep(paste0("\\.", szi, "$"), colnames(padj(cobra)))]
    tth <- truth(cobra)[, grep(paste0("\\.", szi, "$"), colnames(truth(cobra)))]
    colnames(tth) <- gsub(paste0("\\.", szi, "$"), "", colnames(tth))
    
    ## Subset to only the genes that were tested
    tth <- subset(tth, !is.na(tested))
    pvals <- pvals[match(rownames(tth), rownames(pvals)), ]
    padjs <- padjs[match(rownames(tth), rownames(padjs)), ]
    
    ## Set unobserved p-values and adjusted p-values to 1, since these will be
    ## considered "nonsignificant"
    pvals[is.na(pvals)] <- 1
    padjs[is.na(padjs)] <- 1
    
    ## Number of methods calling each variable significant
    tmp1 <- rowSums(padjs <= 0.05)
    tth$nbr_methods <- tmp1[match(rownames(tth), names(tmp1))]
    
    ## Merge adjusted p-values with gene characteristics
    df2 <- dplyr::full_join(melt(as.matrix(padjs)), cbind(Var1 = rownames(tth), tth), 
                            by = "Var1")
    df2$sign <- df2$value <= 0.05
    
    ## Remove extension and data set information from DE method name
    df2$method <- gsub(exts, "", gsub(paste0("\\.", szi, "$"), "", df2$Var2))
    names(colvec) <- gsub(exts, "", names(colvec))
    
    ## Determine which characteristics to logtransform
    allasp <- intersect(
      c("avetpm", "avecount", "vartpm", "fraczero", "fraczeroround",
        "fraczerodiff", "cvtpm", "avecpm", "varcpm", "cvcpm", 
        "fracimputedup", "fracimputeddown", "asinh_nb_minus_zinb_aic"),
      colnames(df2))
    dolog2 <- intersect(c("avetpm", "avecount", "vartpm", "avecpm", "varcpm"),
                        colnames(df2))

    ## For each gene characteristic, populate summary_data with statistics
    ## comparing significant and non-significant genes
    for (y in allasp) {
      message(y)
      nn <- switch(y,
                   avetpm = "average TPM",
                   avecpm = "average CPM",
                   avecount = "average count",
                   vartpm = "variance(TPM)",
                   varcpm = "variance(CPM)",
                   fraczero = "Fraction zeros",
                   fraczeroround = "Fraction zeros after rounding",
                   fraczerodiff = "Difference in fraction zeros between conditions",
                   cvtpm = "coefficient of variation (TPM)",
                   cvcpm = "coefficient of variation (CPM)",
                   fracimputedup = "Fraction of values imputed upwards",
                   fracimputeddown = "Fraction of values imputed downwards", 
                   asinh_nb_minus_zinb_aic = "arcsinh(AIC[NB]-AIC[ZINB])")
      
      ## Log2-transform if necessary
      if (y %in% dolog2) 
        df3 <- df2 %>% dplyr::mutate_(newy = interp(~log2(x), x = as.name(y)))
      else 
        df3 <- df2 %>% dplyr::mutate_(newy = interp(~x, x = as.name(y)))
      
      ## Calculate differences between significant and non-significant genes
      df3 <- df3 %>% dplyr::group_by(Var2) %>% 
        dplyr::summarise(
          tstat = tryCatch(tstatmod(newy[sign], newy[!sign], lmin = 5), 
                           error = function(e) NA),
          snr = tryCatch(snr(newy[sign], newy[!sign], lmin = 5),
                         error = function(e) NA),
          mediandiff = tryCatch(mediandiff(newy[sign], newy[!sign], lmin = 5), 
                                error = function(e) NA))
      
      if (y %in% dolog2) {
        df3 <- df3 %>% dplyr::mutate(charac = paste0("log2_", y))
        pname <- paste0("log2(", y, ")")
      } else { 
        df3 <- df3 %>% dplyr::mutate(charac = y)
        pname <- y
      }
      
      summary_data$stats_charac <- rbind(summary_data$stats_charac, df3)

      ## Plot
      p <- ggplot(df2, aes_string(x = "method", y = pname, 
                                  fill = "method", dodge = "sign", alpha = "sign")) + 
        geom_boxplot() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
        scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        guides(alpha = guide_legend(override.aes = 
                                      list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                           colour = NA))) + 
        stat_summary(fun.data = function(x) {return(c(y = ifelse(y %in% dolog2, 
                                                                 log2(max(df2[, y])),
                                                                 max(df2[, y])),
                                                      label = length(x)))}, 
                     geom = "text", alpha = 1, size = 2, vjust = -1, 
                     position = position_dodge(width = 0.75)) + 
        xlab("") + 
        ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi))
      if (y %in% dolog2) p <- p + ylab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn)))
      else p <- p + ylab(nn)
      print(p)
                                       
      p <- ggplot(df2, aes_string(x = pname, fill = "method", alpha = "sign")) + 
        geom_density() + theme_bw() + scale_fill_manual(values = colvec, name = "") + 
        scale_alpha_manual(values = c(0.2, 0.8), name = "FDR <= 0.05") + 
        facet_wrap(~method, scales = "free_y") + 
        theme(axis.title.y = element_text(size = 13)) + 
        guides(alpha = guide_legend(override.aes = 
                                      list(fill = hcl(c(15, 195), 100, 0, alpha = c(0.2, 0.8)),
                                           colour = NA))) + 
        ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi))
      if (y %in% dolog2) p <- p + xlab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn)))
      else p <- p + xlab(nn)
      print(p)

      p <- ggplot(tth, aes_string(x = "nbr_methods", y = pname, 
                                  group = "nbr_methods")) + 
        geom_boxplot() + theme_bw() + 
        xlab("Number of methods calling gene significant") + 
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 13)) + 
        ggtitle(paste0(ifelse(exts == "", "", paste0(gsub("^_", "", exts), ", ")), szi))
      if (y %in% dolog2) p <- p + ylab(substitute(paste(log[2], "(", nn, ")"), list(nn = nn)))
      else p <- p + ylab(nn)
      print(p)
      
    }
  }
  return(invisible(summary_data))
}

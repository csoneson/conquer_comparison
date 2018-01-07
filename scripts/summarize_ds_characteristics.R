summarize_ds_characteristics <- function(figdir, datasets, exts, dtpext, cols,
                                         singledsfigdir, cobradir, concordancedir, 
                                         dschardir, origvsmockdir, distrdir, 
                                         plotmethods, dstypes, pch_ncells) {
 
  thm <- function() {
    theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13))
  }
  
  ## Read data set characteristics (unfiltered)
  X <- do.call(rbind, lapply(datasets, function(ds) {
    keepgroups <- fromJSON(file = paste0("config/", ds, ".json"))$keepgroups
    readRDS(paste0(dschardir, "/", ds, "_dataset_characteristics_summary_data.rds"))$char_cells_m %>%
      dplyr::filter(ncells == paste0(max(as.numeric(as.character(gsub(" cells per group", "", ncells)))), 
                                     " cells per group")) %>%
      dplyr::filter(mtype %in% c("fraczero", "libsize", "silhouette")) %>%
      dplyr::mutate(condition = condition == keepgroups[1])
  }))
  
  ## Read AIC values
  Y <- do.call(rbind, lapply(datasets, function(ds) {
    do.call(rbind, lapply(exts, function(e) {
      if (file.exists(paste0(distrdir, "/", ds, "mock", e, "_distribution_fit_summary_data.rds"))) {
        readRDS(paste0(distrdir, "/", ds, "mock", e, "_distribution_fit_summary_data.rds"))$GOF_res %>%
          dplyr::select(nbinom_standard_aic, zifnbinom_standard_aic) %>%
          dplyr::mutate(tmp = nbinom_standard_aic - zifnbinom_standard_aic) %>%
          dplyr::mutate(nbinom_minus_zifnbinom_standard_aic = asinh(tmp)) %>%
          dplyr::mutate(dataset = paste0(ds, "mock"), filt = gsub("^_", "", e), 
                        ds = paste0(ds, "mock", e))
      }
    }))
  }))
  Y <- dplyr::left_join(Y, dstypes, by = "dataset") %>%
    dplyr::mutate(ds = gsub("mock", "null", ds))
  
  Y0 <- Y %>% dplyr::filter(!is.na(nbinom_minus_zifnbinom_standard_aic))
  p0 <- ggplot(Y0, aes(x = ds, y = nbinom_minus_zifnbinom_standard_aic)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(alpha = 0.2, position = position_jitter(width = 0.2), size = 0.25) + 
    geom_violin(aes(color = dtype), trim = TRUE, scale = "width", 
                draw_quantiles = 0.5, adjust = 0.5) + 
    guides(fill = FALSE, color = FALSE) + 
    thm() + 
    xlab("") + ylab(expression(arcsinh(AIC[NB] - AIC[ZINB]))) + 
    scale_color_manual(values = c("#B17BA6", "#90C987")) + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = max(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE),
                        label = paste0("n=", sum(!is.na(x)), "\n", 
                                       round(100 * mean(x > 0), 1), "%")))}, 
      geom = "text", alpha = 1, size = 3, vjust = -1) + 
    expand_limits(y = c(min(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE),
                        1.2 * max(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE)))
  
  p0nofilt <- ggplot(Y0 %>% dplyr::filter(filt == ""), 
                     aes(x = ds, y = nbinom_minus_zifnbinom_standard_aic)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(alpha = 0.2, position = position_jitter(width = 0.2), size = 0.25) + 
    geom_violin(aes(color = dtype), trim = TRUE, scale = "width", 
                draw_quantiles = 0.5, adjust = 0.5) + 
    guides(fill = FALSE, color = FALSE) + 
    thm() + 
    xlab("") + ylab(expression(arcsinh(AIC[NB] - AIC[ZINB]))) + 
    scale_color_manual(values = c("#B17BA6", "#90C987")) + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = max(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE),
                        label = paste0("n=", sum(!is.na(x)), "\n", 
                                       round(100 * mean(x > 0), 1), "%")))}, 
      geom = "text", alpha = 1, size = 3, vjust = -1) + 
    expand_limits(y = c(min(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE),
                        1.2 * max(Y0$nbinom_minus_zifnbinom_standard_aic, na.rm = TRUE)))
  
  pdf(paste0(figdir, "/ds_characteristics_final", dtpext, ".pdf"), width = 18, height = 15.6)
  
  X <- X %>% dplyr::filter(!is.na(value))
  X$value[X$mtype == "libsize"] <- log10(X$value[X$mtype == "libsize"])
  p1 <- ggplot(X %>% dplyr::filter(mtype == "fraczero"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + 
    xlab("") + ylab("Fraction zeros per cell") + thm() + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 1.02,
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 3, vjust = -1) + 
    expand_limits(y = 1.06)
  p2 <- ggplot(X %>% dplyr::filter(mtype == "libsize"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + #scale_y_log10() + 
    xlab("") + ylab("log10(library size per cell)") + thm() + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 7.51,
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 3, vjust = -1) + 
    expand_limits(y = 7.75)
  p3 <- ggplot(X %>% dplyr::filter(mtype == "silhouette"), 
               aes(x = dataset, y = value)) + 
    geom_boxplot(outlier.size = -1) +  
    geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(color = condition)) + 
    scale_color_manual(values = structure(c("blue", "red"), names = c(TRUE, FALSE))) + 
    guides(color = FALSE) + 
    xlab("") + ylab("Silhouette width per cell") + thm() + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = 0.39,
                        label = paste0("n=", sum(!is.na(x)))))}, 
      geom = "text", alpha = 1, color = "black", size = 3, vjust = -1) + 
    expand_limits(y = 0.42)

  print(plot_grid(plot_grid(p1, p2, p3, labels = c("A", "B", "C"), align = "h", 
                            rel_widths = c(1, 1, 1), nrow = 1),
                  p0, labels = c("", "D"), rel_heights = c(1, 1.3), ncol = 1))
  
  dev.off()
  
  if (dtpext == "_real") {
    png(paste0(figdir, "/ds_characteristics_for_slides", dtpext, "_1.png"), 
        width = 9, height = 5.6, unit = "in", res = 500)
    print(plot_grid(p2, p1, align = "h", rel_widths = c(1, 1), nrow = 1))
    dev.off()
    
    png(paste0(figdir, "/ds_characteristics_for_slides", dtpext, "_2.png"), 
        width = 9, height = 4.6, unit = "in", res = 500)
    print(p0nofilt)
    dev.off()
  }
}


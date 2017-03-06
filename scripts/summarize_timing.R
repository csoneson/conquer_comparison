summarize_timing <- function(figdir, datasets, exts, dtpext, cols = cols) {
  ## ------------------------------- Timing ----------------------------------- ##
  pdf(paste0(figdir, "/summary_timing", exts, dtpext, ".pdf"), width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/timing/", ds, exts, 
                   "_timing_summary_data.rds"))
  })
  y <- lapply(summary_data_list, function(m) {
    if (!is.null(m$timing_full)) {
      m$timing_full %>% group_by(dataset, filt, ncells, repl) %>%
        dplyr::mutate(timing = timing/max(timing))
    } else {
      NULL
    }
  })
  y <- do.call(rbind, y)
  
  ## Remove extension from method name
  y$method <- gsub(exts, "", y$method)
  y$ncells <- factor(y$ncells, levels = unique(sort(as.numeric(as.character(y$ncells)))))
  
  ## Boxplots
  print(ggplot(y, aes(x = method, y = timing, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = ncells)) + 
          theme_bw() + xlab("") + ylab("Relative timing") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") +
          scale_shape_discrete(name = "Number of cells") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")))

  print(ggplot(y, aes(x = method, y = timing, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2), aes(shape = ncells)) + 
          theme_bw() + xlab("") + ylab("Relative timing") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") +
          scale_shape_discrete(name = "Number of cells") + 
          scale_y_log10() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)) + 
          guides(color = guide_legend(ncol = 2, title = ""),
                 shape = guide_legend(ncol = 2, title = "Number of \ncells per group")))

  ## Barplots
  y %>% group_by(method) %>% dplyr::summarize(mean = mean(timing), sd = sd(timing)) %>%
    ggplot(aes(x = method, y = mean, fill = method)) + 
    geom_errorbar(aes(ymin = min(mean)/2, ymax = mean + sd), width = 0.2) + 
    geom_bar(stat = "identity") + 
    theme_bw() + xlab("") + ylab("Relative timing") + 
    scale_fill_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ## Dependence on number of samples
  y2 <- lapply(summary_data_list, function(m) {
    if (!is.null(m$timing)) {
      m$timing
    } else {
      NULL
    }
  })
  y2 <- do.call(rbind, y2)
  ## Remove extension from method name
  y2$method <- gsub(exts, "", y2$method)
  
  print(y2 %>% group_by(method, dataset, filt) %>%
          arrange(ncells) %>% 
          dplyr::mutate(dt = c(0, diff(timing)),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, diff(ncells)),
                        s = c(0, ncells[1:(length(ncells) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  print(y2 %>% group_by(method, dataset, filt) %>%
          arrange(ncells) %>% 
          dplyr::mutate(dt = c(0, timing[2:(length(timing))]),
                        t = c(0, timing[1:(length(timing) - 1)]),
                        ds = c(0, ncells[2:(length(ncells))]),
                        s = c(0, ncells[1:(length(ncells) - 1)])) %>% 
          dplyr::mutate(reltime = (dt/t)/(ds/s)) %>%
          filter(!is.na(reltime)) %>%
          ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + ylab("Relative change in time per relative increase in number of cells") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  ## Dependence on number of genes
  ngenes <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, 
                   "_ngenes.rds"))
  })
  ngenes <- do.call(rbind, ngenes) %>% group_by(dataset, filt, ncells) %>%
    dplyr::summarize(ngenes = median(ngenes)) %>% 
    mutate(ncells = as.numeric(ncells))
  y3 <- inner_join(y2, ngenes)
  
  ns <- y3 %>% filter(method == y3$method[1]) %>% group_by(ncells) %>% tally()
  ns_keep <- ns$ncells[ns$n > 1]
  
  if (length(ns_keep) > 0) {
    print(y3 %>% filter(ncells %in% ns_keep) %>% group_by(method, ncells, filt) %>%
            dplyr::arrange(ngenes) %>%
            dplyr::filter(length(timing) > 1) %>%
            dplyr::mutate(dt = c(0, timing[2:(length(timing))]),
                          t = c(0, timing[1:(length(timing) - 1)]),
                          dg = c(0, ngenes[2:(length(ngenes))]),
                          g = c(0, ngenes[1:(length(ngenes)) - 1])) %>%
            dplyr::mutate(reltime = (dt/t)/(dg/g)) %>%
            filter(!is.na(reltime)) %>%
            ggplot(aes(x = method, y = reltime, color = method)) + geom_boxplot(outlier.size = -1) + 
            geom_point(position = position_jitter(width = 0.2)) + 
            theme_bw() + xlab("") + 
            ylab("Relative change in time per relative increase in number of genes") + 
            scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_text(size = 13)))
  }
  
  ## 2d heatmap
  ngenes <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/cobra_data/", ds, exts, 
                   "_ngenes.rds"))
  })
  ngenes <- do.call(rbind, ngenes) %>% ungroup() %>% mutate(ncells = as.numeric(ncells)) %>%
    mutate(repl = as.numeric(repl))
  y4 <- lapply(summary_data_list, function(m) {
   m$timing_full
  })
  y4 <- do.call(rbind, y4) %>% mutate(repl = as.numeric(repl))
  y4 <- inner_join(y4, ngenes)
  ## Remove extension from method name
  y4$method <- gsub(exts, "", y4$method)
  
  rbPal <- colorRampPalette(c('red','blue'))
  for (m in unique(y4$method)) {
    pl3d <- scatterplot3d(subset(y4, method == m)$ncells, xlab = "Number of cells per group",  
                          subset(y4, method == m)$ngenes, ylab = "Number of genes", 
                          subset(y4, method == m)$timing, zlab = "Timing", 
                          type = "h", pch = 19, main = m,
                          color = rbPal(20)[as.numeric(cut(subset(y4, method == m)$timing, breaks = 20))])
    model  <- lm(subset(y4, method == m)$timing ~ 
                   subset(y4, method == m)$ncells + subset(y4, method == m)$ngenes)
    pl3d$plane3d(model)
    
    intp <- akima::interp(x = subset(y4, method == m)$ncells, 
                          y = subset(y4, method == m)$ngenes,
                          z = subset(y4, method == m)$timing,
                          duplicate = "mean",
                          xo = seq(min(subset(y4, method == m)$ncells),
                                   max(subset(y4, method == m)$ncells),
                                   length.out = 250),
                          yo = seq(min(subset(y4, method == m)$ngenes),
                                   max(subset(y4, method == m)$ngenes),
                                   length.out = 251))
    rownames(intp$z) <- seq(min(subset(y4, method == m)$ncells),
                            max(subset(y4, method == m)$ncells),
                            length.out = 250)
    colnames(intp$z) <- seq(min(subset(y4, method == m)$ngenes),
                            max(subset(y4, method == m)$ngenes),
                            length.out = 251)
    lab_rows <- rep("", nrow(intp$z))
    w <- sapply(unique(y4$ncells), function(i) which.min(abs(i - as.numeric(rownames(intp$z)))))
    lab_rows[w] <- round(as.numeric(rownames(intp$z))[w])
    lab_cols <- rep("", ncol(intp$z))
    w <- round(seq(1, length(lab_cols), length.out = 5))
    lab_cols[w] <- round(as.numeric(colnames(intp$z))[w])
    pheatmap::pheatmap(intp$z, cluster_rows = FALSE, cluster_cols = FALSE, 
                       labels_row = lab_rows, labels_col = lab_cols,
                       main = m, xlab = "Number of genes")
    
  }
  
  ## Fit power model
  y4$ngenes_cat <- Hmisc::cut2(y4$ngenes, g = 10)
  ## Calculate exponent for single variable (keeping the other fixed)
  calc_expn <- function(y, x) {
    if (length(x) > 3) {
      lm(log(y) ~ log(x))$coef[2]
    } else {
      NA
    }
  }
  ## Fit power model with two predictors
  calc_expn2 <- function(z, x, y) {
    biexp <- function(p, x, y, z) {
      ax <- p[1]
      px <- p[2]
      ay <- p[3]
      py <- p[4]
      sum((z - (ax * x^(px) + ay * y^(py)))^2)
    }
    res <- optim(par = c(1, 1, 1, 1), fn = biexp, x = x, y = y, z = z)
    res$par[c(2, 4)]
  }
  print(y4 %>% group_by(method, ncells, filt) %>% 
          dplyr::summarise(expn = calc_expn(timing, ngenes)) %>%
          ggplot(aes(x = method, y = expn, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + 
          ylab("Exponent in power model of timing vs number of genes") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  print(y4 %>% group_by(method, ngenes_cat, filt) %>% 
          dplyr::summarise(expn = calc_expn(timing, ncells)) %>%
          ggplot(aes(x = method, y = expn, color = method)) + geom_boxplot(outlier.size = -1) + 
          geom_point(position = position_jitter(width = 0.2)) + 
          theme_bw() + xlab("") + 
          ylab("Exponent in power model of timing vs number of cells") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  print(y4 %>% group_by(method, filt) %>% 
          dplyr::summarise(expn = calc_expn2(timing, ncells, ngenes)[1]) %>%
          ggplot(aes(x = method, y = expn, color = method)) + 
          geom_point(size = 3) + 
          theme_bw() + xlab("") + 
          ylab("Exponent in bivariate power model of timing vs number of cells") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  print(y4 %>% group_by(method, filt) %>% 
          dplyr::summarise(expn = calc_expn2(timing, ncells, ngenes)[2]) %>%
          ggplot(aes(x = method, y = expn, color = method)) + 
          geom_point(size = 3) + 
          theme_bw() + xlab("") + 
          ylab("Exponent in bivariate power model of timing vs number of genes") + 
          scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols))), name = "") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 13)))
  
  dev.off()
  
}
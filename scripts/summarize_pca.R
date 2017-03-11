plot_pca <- function(x, title_ext = "", stat) {
  for (scl in c(TRUE)) {
    pca <- prcomp(t(x), scale. = scl)
    annot <- data.frame(id = colnames(x), stringsAsFactors = FALSE) %>%
      tidyr::separate(id, into = c("method", "n_samples", "repl", "dataset", "filt"), 
                      sep = "\\.", remove = FALSE)
    ## Remove extension from method name
    annot$method <- gsub(exts, "", annot$method)
    
    for (cpas in list(c(1, 2), c(1, 3), c(2, 3), c(3, 4))) {
      tmp <- dplyr::full_join(annot, 
                              data.frame(pca$x[, cpas], stringsAsFactors = FALSE) %>% 
                                dplyr::mutate(id = rownames(pca$x)), 
                              by = "id")
      print(ggplot(tmp, 
                   aes_string(x = paste0("PC", cpas[1]), y = paste0("PC", cpas[2]), 
                              color = "method", shape = "dataset")) +
              geom_point(size = 3) + theme_bw() + 
              scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols)))) + 
              ggtitle(paste0(stat, ".scale=", scl, title_ext)) + 
              guides(color = guide_legend(ncol = 2, title = ""),
                     shape = guide_legend(ncol = 2, title = "")))
      print(ggplot(tmp %>%
                     dplyr::group_by(method, dataset) %>% 
                     dplyr::summarize_(PCx = paste0("mean(PC", cpas[1], ")"), 
                                       PCy = paste0("mean(PC", cpas[2], ")")), 
                   aes_string(x = paste0("PCx"), y = paste0("PCy"), 
                              color = "method", shape = "dataset")) +
              geom_point(size = 4) + theme_bw() + 
              xlab(paste0("PC", cpas[1])) + ylab(paste0("PC", cpas[2])) + 
              scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols)))) + 
              ggtitle(paste0(stat, ".scale=", scl, title_ext)) + 
              guides(color = guide_legend(ncol = 2, title = ""),
                     shape = guide_legend(ncol = 2, title = "")))
      print(ggplot(tmp %>%
                     dplyr::group_by(method, dataset) %>% 
                     dplyr::summarize_(PCx = paste0("mean(PC", cpas[1], ")"), 
                                       PCy = paste0("mean(PC", cpas[2], ")")), 
                   aes_string(x = paste0("PCx"), y = paste0("PCy"), 
                              color = "method", shape = "dataset",
                              label = "method")) +
              geom_point(size = 4) + theme_bw() + geom_text_repel() + 
              xlab(paste0("PC", cpas[1])) + ylab(paste0("PC", cpas[2])) + 
              scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols)))) + 
              ggtitle(paste0(stat, ".scale=", scl, title_ext)) + 
              guides(color = guide_legend(ncol = 2, title = ""),
                     shape = guide_legend(ncol = 2, title = "")))
      print(ggplot(data.frame(id = rownames(pca$rotation), pca$rotation[, cpas]), 
                   aes_string(x = paste0("PC", cpas[1]), y = paste0("PC", cpas[2]), 
                              label = "id")) + geom_point() + geom_text_repel() +
              geom_segment(aes_string(x = 0, y = 0, 
                                      xend = paste0("PC", cpas[1]), yend = paste0("PC", cpas[2])), 
                           arrow = arrow(length = unit(0.03, "npc")), linetype = "dashed") + 
              theme_bw() + 
              ggtitle(paste0(stat, ".scale=", scl, title_ext)))
      print(ggbiplot(pca, scale = 0, groups = annot$method, ellipse = TRUE,
                     ellipse.prob = 0.68, alpha = 0, var.axes = TRUE, choices = cpas) + 
              theme_bw() + 
              scale_color_manual(values = structure(cols, names = gsub(exts, "", names(cols)))) + 
              ggtitle(paste0(stat, ".scale=", scl, title_ext)) + 
              geom_point(data = tmp %>%
                           dplyr::group_by(method) %>% 
                           dplyr::summarize_(PCx = paste0("mean(PC", cpas[1], ")"), 
                                             PCy = paste0("mean(PC", cpas[2], ")")),
                         aes(x = PCx, y = PCy, col = method)) + 
              guides(color = guide_legend(ncol = 2, title = ""),
                     shape = guide_legend(ncol = 2, title = "")) + 
              coord_equal(ratio = diff(range(pca$x[, cpas[1]]))/diff(range(pca$x[, cpas[2]]))))
    }
  }
}

summarize_pca <- function(figdir, datasets, exts, dtpext, cols = cols) {
  ## ---------------------------------- PCA ----------------------------------- ##
  pdf(paste0(figdir, "/summary_pca", exts, dtpext, ".pdf"), width = 10, height = 7)
  summary_data_list <- lapply(datasets, function(ds) {
    readRDS(paste0("figures/results_characterization/", ds, exts, 
                   "_results_characterization_summary_data.rds"))
  })
  ## PCA of significant gene characteristics
  lapply(c(list(datasets), as.list(datasets)), function(ds) {
    for (stat in c("tstat")) {
      x <- lapply(summary_data_list[ds], function(m) {
        m$stats_charac %>% dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
          dplyr::filter_(paste0("is.finite(", stat, ")")) %>% 
          dplyr::mutate(Var2 = paste0(Var2, ".", dataset, ".", filt)) %>%
          dplyr::select_("Var2", stat, "charac")  %>%
          dplyr::filter(charac != "fraczerodiff") %>%
          dplyr::filter(charac != "log2_avecount")
      })
      x <- do.call(rbind, x) %>% dcast(charac ~ Var2, value.var = stat)
      rownames(x) <- x$charac
      x$charac <- NULL
      x[is.na(x)] <- 0
      
      plot_pca(x, title_ext = "", stat = stat)
    }
  })
  
  ## The same as above, but centering and scaling each statistic within each data set before merging
  lapply(c(list(datasets)), function(ds) {
    for (stat in c("tstat")) {
      x <- lapply(summary_data_list[ds], function(m) {
        m$stats_charac %>% dplyr::filter_(paste0("!is.na(", stat, ")")) %>% 
          dplyr::filter_(paste0("is.finite(", stat, ")")) %>% 
          dplyr::mutate(Var2 = paste0(Var2, ".", dataset, ".", filt)) %>%
          dplyr::select_("Var2", stat, "charac")  %>%
          dplyr::filter(charac != "fraczerodiff") %>%
          dplyr::filter(charac != "log2_avecount") %>%
          dplyr::group_by(charac) %>%
          dplyr::mutate_(newy = interp(~scale(x, scale = FALSE), x = as.name(stat))) %>%
          dplyr::select(Var2, newy, charac)
      })
      x <- do.call(rbind, x) %>% dcast(charac ~ Var2, value.var = "newy")
      rownames(x) <- x$charac
      x$charac <- NULL
      x[is.na(x)] <- 0
      
      plot_pca(x, title_ext = ", after scaling within each data set", stat = stat)
     }
  })
  dev.off()
}
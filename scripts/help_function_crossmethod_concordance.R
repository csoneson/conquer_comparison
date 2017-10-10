suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

## Get all subclusters from an hclust object
get_subclusters <- function(hcl) {
  m <- hcl$merge
  labs <- hcl$labels
  L <- list()
  for (i in seq_len(nrow(m))) {
    tmp <- c()
    if (m[i, 1] < 0) tmp <- c(tmp, labs[-m[i, 1]])
    else tmp <- c(tmp, L[[m[i, 1]]])
    if (m[i, 2] < 0) tmp <- c(tmp, labs[-m[i, 2]])
    else tmp <- c(tmp, L[[m[i, 2]]])
    L[[i]] <- sort(tmp)
  }
  L
}

## Visualize cross-method consistency (average pairwise AUCC between methods for
## top-k0 genes)
help_function_crossmethod_concordance <- function(concordance_betweenmethods_pairwise, 
                                                  k0, titleext = "") {
  plots <- list()
  
  ## Make sure that both combinations of each method pair are represented
  cm1 <- concordance_betweenmethods_pairwise %>% dplyr::ungroup()
  cm2 <- concordance_betweenmethods_pairwise %>% dplyr::ungroup()
  cm2[, c("method1", "method2")] <- cm2[, c("method2", "method1")]
  cm <- rbind(cm1, cm2)
  
  ## Calculate average area under concordance curve across all data set
  ## instances, for each pair of methods
  cmcons <- cm  %>% 
    dplyr::filter(k == k0) %>% dplyr::group_by(method1, method2) %>%
    dplyr::summarize(meanAUCs = mean(AUCs)) %>% as.data.frame()
  ## Add concordance for each method with itself
  cmcons <- rbind(cmcons, data.frame(method1 = unique(c(cmcons$method1, cmcons$method2)),
                                     method2 = unique(c(cmcons$method1, cmcons$method2)),
                                     meanAUCs = 1, stringsAsFactors = FALSE))
  cmcons <- dcast(cmcons, method1 ~ method2, value.var = "meanAUCs")
  rownames(cmcons) <- cmcons$method1
  cmcons$method1 <- NULL
  stopifnot(all(rownames(cmcons) == colnames(cmcons)))
  stopifnot(all((cmcons == t(cmcons))[!is.na(cmcons == t(cmcons))]))
  
  ## Hierarchical clustering based on 1 - cmcons
  tmpdist <- 1 - cmcons
  tmpdist[is.na(tmpdist)] <- 1
  hcl_average <- hclust(as.dist(tmpdist))
  plot(hcl_average)
  ## Get all subclusters
  subclusters_average <- get_subclusters(hcl_average)
  
  ## Get all subclusters for individual data set instances
  cmtmp <- cm %>% dplyr::filter(k == k0) %>% dplyr::mutate(grp = paste(dataset, filt, ncells, repl, sep = "."))
  uniqvals <- unique(cmtmp$grp)
  subclusters_all <- lapply(uniqvals, function(i) {
    cmtmp2 <- cmtmp %>% dplyr::filter(grp == i) %>% dplyr::select(method1, method2, AUCs)
    cmtmp2 <- rbind(cmtmp2, data.frame(method1 = unique(c(cmtmp2$method1, cmtmp2$method2)),
                                       method2 = unique(c(cmtmp2$method1, cmtmp2$method2)),
                                       AUCs = 1, stringsAsFactors = FALSE))
    cmtmp2 <- dcast(cmtmp2, method1 ~ method2, value.var = "AUCs")
    rownames(cmtmp2) <- cmtmp2$method1
    cmtmp2$method1 <- NULL
    stopifnot(all(rownames(cmtmp2) == colnames(cmtmp2)))
    stopifnot(all((cmtmp2 == t(cmtmp2))[!is.na(cmtmp2 == t(cmtmp2))]))
    get_subclusters(hclust(as.dist(1 - cmtmp2)))
  })
  
  ## Get stability values for ech subcluster in subclusters_average
  stability_scores <- rowMeans(sapply(subclusters_all, function(w) {
    subclusters_average %in% w
  }))
  
  plots[["stability_scores"]] <- stability_scores
  plots[["subclusters_average"]] <- subclusters_average
  plots[["clustering_average"]] <- hcl_average
  
  phm <- pheatmap(cmcons, clustering_distance_rows = as.dist(tmpdist), 
                  clustering_distance_cols = as.dist(tmpdist), clustering_method = "complete", 
                  main = paste0("Area under method/method concordance curve,", 
                                "\naveraged across all data set instances, top-", 
                                k0, " genes, ", titleext), 
                  color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                            "RdYlBu")))(100),
                  breaks = seq(0, 1, length.out = 101))

  ## Visualize full distributions of cross-method consistencies across data set
  ## instances
  cmdist <- cm %>% dplyr::filter(k == k0) %>% as.data.frame()
  
  ## Add "diagonal" similarities between a method and itself
  allm <- unique(c(cmdist$method1, cmdist$method2))
  tmp <- cmdist %>% dplyr::group_by(dataset, filt, ncells, repl) %>%
    dplyr::filter(row_number() == 1) %>% ungroup()
  tmp <- tmp[rep(seq_len(nrow(tmp)), length(allm)), ] %>%
    dplyr::mutate(method1 = rep(allm, each = nrow(tmp)), 
                  method2 = rep(allm, each = nrow(tmp)), 
                  AUCs = 1, AUC = k^2/2) %>% as.data.frame()
  cmdist <- rbind(cmdist, tmp)
  
  cmdist <- cmdist %>% dplyr::filter(method1 <= method2)
  
  ## Order methods by hierarchical clustering result and visualize distributions
  ## of AUCs in color code
  hclord <- rev(hcl_average$labels[hcl_average$order])
  ## Order AUCs in increasing order for each method pair
  cmdist2 <- cmdist %>% dplyr::group_by(method1, method2) %>%
    dplyr::arrange(AUCs) %>% dplyr::mutate(ordr = 1:length(AUCs), ycoord = 1)
  for (i in 1:nrow(cmdist2)) {
    if (match(cmdist2[i, "method1"], hclord) > match(cmdist2[i, "method2"], hclord)) {
      cmdist2[i, c("method1", "method2")] <- cmdist2[i, c("method2", "method1")]
    }
  }
  cmdist2$method1 <- factor(cmdist2$method1, levels = hclord[hclord %in% cmdist2$method1])
  cmdist2$method2 <- factor(cmdist2$method2, levels = hclord[hclord %in% cmdist2$method2])

  p <- ggplot(cmdist2, aes(x = ordr, y = ycoord)) + geom_raster(aes(fill = AUCs)) +
    facet_grid(method1 ~ method2, scales = "free_x") + theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          strip.text.x = element_text(size = 12, angle = 90),
          strip.text.y = element_text(size = 12, angle = 0), 
          panel.grid = element_blank(),
          panel.border = element_blank(), 
          strip.background = element_rect(colour = "white"),
          legend.position = c(0.1, 0.3),
          panel.spacing = unit(0.15, "lines")) +
    xlab("") + ylab("") + ggtitle(titleext) + 
    scale_fill_continuous(low = "black", high = "yellow", 
                          name = paste0("AUCC,\ntop-", k0, "\ngenes"), 
                          limits = c(0, 1))
  plots[["concordancedistr_color"]] <- p
  print(p)
  
  ## Plot dependence of AUC on number of cells per group, for each pair of
  ## methods
  cmdist2 <- dplyr::mutate(cmdist2, ncells = as.numeric(as.character(ncells)))
  p <- ggplot(cmdist2, aes(x = ncells, y = AUCs)) + geom_point(size = 0.5) +
    geom_smooth() + ggtitle(titleext) + 
    facet_grid(method1 ~ method2) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
          strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank(),
          strip.background = element_rect(colour = "white"),
          panel.spacing = unit(0.1, "lines")) +
    xlab("Number of cells per group") + ylab(paste0("AUCC, top-", k0, " genes"))
  plots[["concordance_dep_ncells"]] <- p
  print(p)
  
  ## Plot Spearman correlation between number of cells per group and AUCs, for
  ## each pair of methods
  p <- ggplot(cmdist2 %>% dplyr::group_by(method1, method2) %>%
                dplyr::summarize(spearman = cor(ncells, AUCs, method = "spearman")) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(method1 = factor(method1, levels = hclord[hclord %in% method1]),
                              method2 = factor(method2, levels = hclord[hclord %in% method2])),
              aes(x = 1, y = 1)) + geom_raster(aes(fill = spearman)) +
    facet_grid(method1 ~ method2, scales = "free_x") + theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          strip.text.x = element_text(size = 12, angle = 90),
          strip.text.y = element_text(size = 12, angle = 0), 
          panel.grid = element_blank(),
          panel.border = element_blank(), 
          strip.background = element_rect(colour = "white"),
          legend.position = c(0.1, 0.3),
          panel.spacing = unit(0.15, "lines")) +
    xlab("") + ylab("") + ggtitle(titleext) + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                          name = paste0("Spearman\ncorrelation\nbetween\nnumber of cells\nand AUCC,\ntop-",
                                        k0, " genes"), 
                          limits = c(-1, 1))
  plots[["concordance_dep_ncells_color"]] <- p
  print(p)
  
  ## Plot dependence of AUC on average silhouette width, for each pair of
  ## methods
  if ("silhouette_avg" %in% colnames(cmdist2)) {
    p <- ggplot(cmdist2, aes(x = silhouette_avg, y = AUCs)) + geom_point(size = 0.5) +
      geom_smooth() + ggtitle(titleext) + 
      facet_grid(method1 ~ method2) + theme_bw() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank(), strip.text.x = element_text(size = 12, angle = 90),
            strip.text.y = element_text(size = 12, angle = 0), panel.grid = element_blank(),
            strip.background = element_rect(colour = "white"),
            panel.spacing = unit(0.1, "lines")) +
      xlab("Average silhouette width") + ylab(paste0("AUCC, top-", k0, " genes"))
    plots[["concordance_dep_silhouette"]] <- p
    print(p)
    
    ## Plot Spearman correlation between silhouette width and AUCs, for
    ## each pair of methods
    p <- ggplot(cmdist2 %>% dplyr::group_by(method1, method2) %>%
                  dplyr::summarize(spearman = cor(silhouette_avg, AUCs, method = "spearman")) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(method1 = factor(method1, levels = hclord[hclord %in% method1]),
                                method2 = factor(method2, levels = hclord[hclord %in% method2])),
                aes(x = 1, y = 1)) + geom_raster(aes(fill = spearman)) +
      facet_grid(method1 ~ method2, scales = "free_x") + theme_bw() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank(), 
            strip.text.x = element_text(size = 12, angle = 90),
            strip.text.y = element_text(size = 12, angle = 0), 
            panel.grid = element_blank(),
            panel.border = element_blank(), 
            strip.background = element_rect(colour = "white"),
            legend.position = c(0.1, 0.3),
            panel.spacing = unit(0.15, "lines")) +
      xlab("") + ylab("") + ggtitle(titleext) + 
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           name = paste0("Spearman\ncorrelation\nbetween\naverage\nsilhouette ", 
                                         "width\nand AUCC,\ntop-", k0, " genes"), 
                           limits = c(-1, 1))
    plots[["concordance_dep_silhouette_color"]] <- p
    print(p)
  }
  
  invisible(plots)
}
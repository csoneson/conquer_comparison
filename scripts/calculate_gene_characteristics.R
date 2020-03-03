suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(edgeR))

calculate_gene_characteristics <- function(L, do.plot = FALSE, title.ext = "") {
  ## Census counts
  cds <- newCellDataSet(L$tpm, 
                        phenoData = new("AnnotatedDataFrame", 
                                        data = data.frame(condition = L$condt, 
                                                          row.names = colnames(L$tpm))),
                        expressionFamily = tobit())
  censuscounts <- relative2abs(cds)
  
  ## Plot
  if (do.plot) {
    ## Census count distributions
    print(reshape2::melt(censuscounts + 1) %>% dplyr::mutate(condition = L$condt[Var2]) %>%
            ggplot(aes(x = value, group = Var2, color = condition)) + 
            geom_density() + scale_x_log10() + theme_bw() + 
            xlab("Census count + 1") + ggtitle(paste0("Census count distribution per cell", title.ext)))
    print(reshape2::melt(censuscounts) %>% dplyr::mutate(condition = L$condt[Var2]) %>%
            ggplot(aes(x = value, group = Var2, color = condition)) + 
            geom_density() + scale_x_log10() + theme_bw() + 
            xlab("Census count, zeros excluded") + 
            ggtitle(paste0("Census count distribution per cell", title.ext)))
    
    ## Count distributions
    print(reshape2::melt(L$count + 1) %>% dplyr::mutate(condition = L$condt[Var2]) %>% 
            ggplot(aes(x = value, group = Var2, color = condition)) + 
            geom_density() + scale_x_log10() + theme_bw() + 
            xlab("Count + 1") + ggtitle(paste0("Count distribution per cell", title.ext)))
    print(reshape2::melt(L$count) %>% dplyr::mutate(condition = L$condt[Var2]) %>% 
            ggplot(aes(x = value, group = Var2, color = condition)) + 
            geom_density() + scale_x_log10() + theme_bw() + 
            xlab("Count, zeros excluded") + ggtitle(paste0("Count distribution per cell", title.ext)))
  }
  
  avecensuscount <- data.frame(avecensuscount = apply(censuscounts, 1, mean), 
                               gene = rownames(censuscounts))
  fraczerocensus <- data.frame(fraczerocensus = apply(censuscounts, 1, function(x) mean(x == 0)),
                               gene = rownames(censuscounts))
  
  ## CPMs
  cpms <- edgeR::cpm(L$count, lib.size = colSums(L$count) * edgeR::calcNormFactors(L$count))
  avecpm <- data.frame(avecpm = apply(cpms, 1, mean), gene = rownames(cpms))
  varcpm <- data.frame(varcpm = apply(cpms, 1, var), gene = rownames(cpms))
  cvcpm <- data.frame(cvcpm = apply(cpms, 1, sd)/apply(cpms, 1, mean), gene = rownames(cpms))
  
  ## Average raw count
  avecount <- data.frame(avecount = apply(L$count, 1, mean), gene = rownames(L$count))
  
  ## Average TPM
  avetpm <- data.frame(avetpm = apply(L$tpm, 1, mean), gene = rownames(L$tpm))
  
  ## Fraction zeros
  fraczero <- data.frame(fraczero = apply(L$count, 1, function(x) mean(x == 0)),
                         fraczeroround = apply(round(L$count), 1, function(x) mean(x == 0)), 
                         fraczero1 = apply(L$count[, L$condt == levels(factor(L$condt))[1]], 
                                           1, function(x) mean(x == 0)),
                         fraczero2 = apply(L$count[, L$condt == levels(factor(L$condt))[2]], 
                                           1, function(x) mean(x == 0)), 
                         gene = rownames(L$count))
  fraczero$fraczerodiff <- abs(fraczero$fraczero1 - fraczero$fraczero2)
  
  ## Variance of TPMs
  vartpm <- data.frame(vartpm = apply(L$tpm, 1, var), gene = rownames(L$tpm))
  
  ## Coefficient of variation of TPMs
  cvtpm <- data.frame(cvtpm = apply(L$tpm, 1, sd)/apply(L$tpm, 1, mean), gene = rownames(L$tpm))
  
  ## Put all together
  df2 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
                list(vartpm, fraczero, avecount, avetpm, cvtpm, avecensuscount, 
                     fraczerocensus, avecpm, varcpm, cvcpm))
  
  ## Number of imputed values
  if (!is.null(L$nimp)) {
    n <- ncol(L$count)
    nimp <- L$nimp %>% dplyr::mutate(fracimputedup = nbr_increased/n,
                                     fracimputeddown = nbr_decreased/n,
                                     fracnonimputed = nbr_unchanged/n) %>%
      dplyr::select(gene, fracimputedup, fracimputeddown, fracnonimputed)
    df2 <- merge(df2, nimp, by = "gene", all = TRUE)
  } else {
    df2$fracimputedup <- 0
    df2$fracimputeddown <- 0
    df2$fracnonimputed <- 1
  }
  
  ## Add column "tested", which is TRUE for all genes (all genes are sent into the test)
  df2$tested <- TRUE

  list(characs = df2)
}
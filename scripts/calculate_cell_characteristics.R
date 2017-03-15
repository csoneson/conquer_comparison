suppressPackageStartupMessages(library(monocle))

calculate_cell_characteristics <- function(L) {
  cds <- newCellDataSet(L$tpm, 
                        phenoData = new("AnnotatedDataFrame", 
                                        data = data.frame(condition = L$condt, 
                                                          row.names = colnames(L$tpm))))
  censuscounts <- relative2abs(cds)
  
  libsize <- data.frame(libsize = colSums(L$count), cell = colnames(L$count))
  libsizecensus <- data.frame(libsizecensus = colSums(censuscounts), 
                              cell = colnames(censuscounts))
  fraczero <- data.frame(fraczero = colMeans(L$count == 0), cell = colnames(L$count))
  fraczeroround <- data.frame(fraczeroround = colMeans(round(L$count) == 0), 
                              cell = colnames(L$count))
  fraczerocensus <- data.frame(fraczerocensus = colMeans(censuscounts == 0), 
                               cell = colnames(censuscounts))
  df3 <- Reduce(function(...) dplyr::full_join(..., by = "cell"),
                list(libsize, fraczero, fraczeroround, libsizecensus, fraczerocensus))
  
  list(characs = df3)
}
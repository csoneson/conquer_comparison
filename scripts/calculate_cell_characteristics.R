suppressPackageStartupMessages(library(monocle))

calculate_cell_characteristics <- function(L) {
  cds <- newCellDataSet(L$tpm, 
                        phenoData = new("AnnotatedDataFrame", 
                                        data = data.frame(condition = L$condt, 
                                                          row.names = colnames(L$tpm))))
  censuscounts <- relative2abs(cds)
  
  libsize <- data.frame(libsize = colSums(L$count), cell = colnames(L$count))
  libsizecensus <- data.frame(libsizecensus = colSums(censuscounts), cell = colnames(censuscounts))
  fraczerocell <- data.frame(fraczero = colMeans(L$count == 0), cell = colnames(L$count))
  fraczerocellcensus <- data.frame(fraczerocensus = colMeans(censuscounts == 0), 
                                   cell = colnames(censuscounts))
  df3 <- Reduce(function(...) merge(..., by = "cell", all = TRUE),
                list(libsize, fraczerocell, libsizecensus, fraczerocellcensus))
  
  list(characs = df3)
}
suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596imputemock.rds", 
          subfile = "subsets/GSE74596imputemock_subsets.rds",
          resfilebase = "results/GSE74596imputemock",
          figfilebase = "figures/diffexpression/GSE74596imputemock", 
          groupid = "source_name_ch1", 
          keepgroups = "Single_cell_RNA-seq_NKT0", 
          seed = 42, 
          sizes = c(22, 12, 6), 
          nreps = c(1, 5, 5),
          impute = "yes")
write(toJSON(L), file = "config/GSE74596imputemock.json")

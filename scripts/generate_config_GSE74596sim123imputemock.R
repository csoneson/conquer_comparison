suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596sim123imputemock.rds", 
          subfile = "subsets/GSE74596sim123imputemock_subsets.rds",
          resfilebase = "results/GSE74596sim123imputemock",
          figfilebase = "figures/diffexpression/GSE74596sim123imputemock", 
          groupid = "group", 
          keepgroups = "Single_cell_RNA-seq_NKT0", 
          seed = 42, 
          sizes = c(22, 12, 6), 
          nreps = c(1, 5, 5), 
          impute = "yes")
write(toJSON(L), file = "config/GSE74596sim123imputemock.json")

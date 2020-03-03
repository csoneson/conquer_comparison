suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596sim123scimputemock.rds", 
          subfile = "subsets/GSE74596sim123scimputemock_subsets.rds",
          resfilebase = "results/GSE74596sim123scimputemock",
          figfilebase = "figures/diffexpression/GSE74596sim123scimputemock", 
          groupid = "group", 
          keepgroups = "Single_cell_RNA-seq_NKT0", 
          seed = 42, 
          sizes = c(22, 12, 6), 
          nreps = c(1, 5, 5), 
          impute = "scimpute")
write(toJSON(L), file = "config/GSE74596sim123scimputemock.json")

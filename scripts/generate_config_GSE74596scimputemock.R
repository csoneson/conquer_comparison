suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596scimputemock.rds", 
          subfile = "subsets/GSE74596scimputemock_subsets.rds",
          resfilebase = "results/GSE74596scimputemock",
          figfilebase = "figures/diffexpression/GSE74596scimputemock", 
          groupid = "source_name_ch1", 
          keepgroups = "Single_cell_RNA-seq_NKT0", 
          seed = 42, 
          sizes = c(22, 12, 6), 
          nreps = c(1, 5, 5),
          impute = "scimpute")
write(toJSON(L), file = "config/GSE74596scimputemock.json")

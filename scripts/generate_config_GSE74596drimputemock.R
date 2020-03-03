suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596drimputemock.rds", 
          subfile = "subsets/GSE74596drimputemock_subsets.rds",
          resfilebase = "results/GSE74596drimputemock",
          figfilebase = "figures/diffexpression/GSE74596drimputemock", 
          groupid = "source_name_ch1", 
          keepgroups = "Single_cell_RNA-seq_NKT0", 
          seed = 42, 
          sizes = c(22, 12, 6), 
          nreps = c(1, 5, 5),
          impute = "drimpute")
write(toJSON(L), file = "config/GSE74596drimputemock.json")

suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596drimpute.rds", 
          subfile = "subsets/GSE74596drimpute_subsets.rds",
          resfilebase = "results/GSE74596drimpute",
          figfilebase = "figures/diffexpression/GSE74596drimpute", 
          groupid = "source_name_ch1", 
          keepgroups = c("Single_cell_RNA-seq_NKT0",
                         "Single_cell_RNA-seq_NKT17"), 
          seed = 42, 
          sizes = c(44, 22, 12, 6), 
          nreps = c(1, 5, 5, 5),
          impute = "drimpute")
write(toJSON(L), file = "config/GSE74596drimpute.json")

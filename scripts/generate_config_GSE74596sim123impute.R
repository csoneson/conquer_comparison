suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596sim123impute.rds", 
          subfile = "subsets/GSE74596sim123impute_subsets.rds",
          resfilebase = "results/GSE74596sim123impute",
          figfilebase = "figures/diffexpression/GSE74596sim123impute", 
          groupid = "group", 
          keepgroups = c("Single_cell_RNA-seq_NKT0",
                         "Single_cell_RNA-seq_NKT17"), 
          seed = 42, 
          sizes = c(44, 22, 12, 6), 
          nreps = c(1, 5, 5, 5),
          impute = "yes")
write(toJSON(L), file = "config/GSE74596sim123impute.json")

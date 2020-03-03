suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596sim123drimpute.rds", 
          subfile = "subsets/GSE74596sim123drimpute_subsets.rds",
          resfilebase = "results/GSE74596sim123drimpute",
          figfilebase = "figures/diffexpression/GSE74596sim123drimpute", 
          groupid = "group", 
          keepgroups = c("Single_cell_RNA-seq_NKT0",
                         "Single_cell_RNA-seq_NKT17"), 
          seed = 42, 
          sizes = c(44, 22, 12, 6), 
          nreps = c(1, 5, 5, 5),
          impute = "drimpute")
write(toJSON(L), file = "config/GSE74596sim123drimpute.json")

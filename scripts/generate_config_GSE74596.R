suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE74596.rds", 
          subfile = "subsets/GSE74596_subsets.rds",
          resfilebase = "results/GSE74596",
          figfilebase = "figures/diffexpression/GSE74596", 
          groupid = "source_name_ch1", 
          keepgroups = c("Single_cell_RNA-seq_NKT0",
                         "Single_cell_RNA-seq_NKT17"), 
          seed = 42, 
          sizes = c(44, 22, 12, 6), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE74596.json")

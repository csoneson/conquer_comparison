suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/EMTAB2805.rds", 
          subfile = "subsets/EMTAB2805_subsets.rds",
          resfilebase = "results/EMTAB2805",
          figfilebase = "figures/diffexpression/EMTAB2805", 
          groupid = "cell_cycle_stage", 
          keepgroups = NULL, 
          seed = 42, 
          sizes = c(96, 48, 24, 12), 
          nreps = c(1, 3, 3, 3))
write(toJSON(L), file = "config/EMTAB2805.json")

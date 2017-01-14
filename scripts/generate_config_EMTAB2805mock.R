suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/EMTAB2805mock.rds", 
          subfile = "subsets/EMTAB2805mock_subsets.rds",
          resfilebase = "results/EMTAB2805mock",
          figfilebase = "figures/diffexpression/EMTAB2805mock", 
          groupid = "cell_cycle_stage", 
          keepgroups = "G2M", 
          seed = 42, 
          sizes = c(48, 24, 12, 6), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/EMTAB2805mock.json")

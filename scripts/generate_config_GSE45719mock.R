suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE45719.rds", 
          subfile = "subsets/GSE45719mock_subsets.rds",
          resfilebase = "results/GSE45719mock",
          figfilebase = "figures/diffexpression/GSE45719mock", 
          groupid = "source_name_ch1", 
          keepgroups = c("16-cell stage blastomere"), 
          seed = 42, 
          sizes = c(24, 14, 12, 6), 
          nreps = c(1, 3, 3, 3))
write(toJSON(L), file = "config/GSE45719mock.json")

suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE60749-GPL13112mock.rds", 
          subfile = "subsets/GSE60749-GPL13112mock_subsets.rds",
          resfilebase = "results/GSE60749-GPL13112mock",
          figfilebase = "figures/diffexpression/GSE60749-GPL13112mock", 
          groupid = c("source_name_ch1", "characteristics_ch1.1"), 
          keepgroups = "v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF", 
          seed = 42, 
          sizes = c(47, 24, 12, 6), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE60749-GPL13112mock.json")

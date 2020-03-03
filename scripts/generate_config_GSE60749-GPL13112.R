suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE60749-GPL13112.rds", 
          subfile = "subsets/GSE60749-GPL13112_subsets.rds",
          resfilebase = "results/GSE60749-GPL13112",
          figfilebase = "figures/diffexpression/GSE60749-GPL13112", 
          groupid = c("source_name_ch1", "characteristics_ch1.1"), 
          keepgroups = c("v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF",
                         "v6.5 mouse embryonic stem cells.culture conditions: serum+LIF"), 
          seed = 42, 
          sizes = c(90, 48, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE60749-GPL13112.json")

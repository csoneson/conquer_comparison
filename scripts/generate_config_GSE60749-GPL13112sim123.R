suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE60749-GPL13112sim123.rds", 
          subfile = "subsets/GSE60749-GPL13112sim123_subsets.rds",
          resfilebase = "results/GSE60749-GPL13112sim123",
          figfilebase = "figures/diffexpression/GSE60749-GPL13112sim123", 
          groupid = "group", 
          keepgroups = c("v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF",
                         "v6.5 mouse embryonic stem cells.culture conditions: serum+LIF"), 
          seed = 42, 
          sizes = c(90, 48, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE60749-GPL13112sim123.json")

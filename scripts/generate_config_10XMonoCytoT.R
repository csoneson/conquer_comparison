suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/10XMonoCytoT.rds", 
          subfile = "subsets/10XMonoCytoT_subsets.rds",
          resfilebase = "results/10XMonoCytoT",
          figfilebase = "figures/diffexpression/10XMonoCytoT", 
          groupid = "group", 
          keepgroups = c("monocyte", "cytotoxict"), 
          seed = 42, 
          sizes = c(240, 120, 48, 24), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/10XMonoCytoT.json")

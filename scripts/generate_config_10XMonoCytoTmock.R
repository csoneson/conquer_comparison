suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/10XMonoCytoTmock.rds", 
          subfile = "subsets/10XMonoCytoTmock_subsets.rds",
          resfilebase = "results/10XMonoCytoTmock",
          figfilebase = "figures/diffexpression/10XMonoCytoTmock", 
          groupid = "group", 
          keepgroups = c("cytotoxict"), 
          seed = 42, 
          sizes = c(120, 48, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/10XMonoCytoTmock.json")

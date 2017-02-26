suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/EGEUV1mock.rds", 
          subfile = "subsets/EGEUV1mock_subsets.rds",
          resfilebase = "results/EGEUV1mock",
          figfilebase = "figures/diffexpression/EGEUV1mock", 
          groupid = c("center", "population"), 
          keepgroups = c("University of Geneva.CEU"), 
          seed = 42, 
          sizes = c(44, 24, 12, 6), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/EGEUV1mock.json")

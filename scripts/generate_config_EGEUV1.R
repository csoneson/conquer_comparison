suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/EGEUV1.rds", 
          subfile = "subsets/EGEUV1_subsets.rds",
          resfilebase = "results/EGEUV1",
          figfilebase = "figures/diffexpression/EGEUV1", 
          groupid = c("center", "population"), 
          keepgroups = c("University of Geneva.CEU",
                         "University of Geneva.YRI"), 
          seed = 42, 
          sizes = c(87, 48, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/EGEUV1.json")

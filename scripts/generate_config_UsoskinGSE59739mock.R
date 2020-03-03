suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/UsoskinGSE59739mock.rds", 
          subfile = "subsets/UsoskinGSE59739mock_subsets.rds",
          resfilebase = "results/UsoskinGSE59739mock",
          figfilebase = "figures/diffexpression/UsoskinGSE59739mock", 
          groupid = "Picking.sessions.Level.3", 
          keepgroups = "RT-1.NP1", 
          seed = 42, 
          sizes = c(24, 12, 6), 
          nreps = c(1, 5, 5))
write(toJSON(L), file = "config/UsoskinGSE59739mock.json")

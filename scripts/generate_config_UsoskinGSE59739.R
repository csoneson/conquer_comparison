suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/UsoskinGSE59739.rds", 
          subfile = "subsets/UsoskinGSE59739_subsets.rds",
          resfilebase = "results/UsoskinGSE59739",
          figfilebase = "figures/diffexpression/UsoskinGSE59739", 
          groupid = "Picking.sessions.Level.3", 
          keepgroups = c("RT-1.NP1", "RT-1.TH"), 
          seed = 42, 
          sizes = c(58, 36, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/UsoskinGSE59739.json")

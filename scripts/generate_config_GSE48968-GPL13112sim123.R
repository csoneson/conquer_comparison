suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE48968-GPL13112sim123.rds", 
          subfile = "subsets/GSE48968-GPL13112sim123_subsets.rds",
          resfilebase = "results/GSE48968-GPL13112sim123",
          figfilebase = "figures/diffexpression/GSE48968-GPL13112sim123", 
          groupid = "group", 
          keepgroups = c("BMDC (1h LPS Stimulation)",
                         "BMDC (4h LPS Stimulation)"), 
          seed = 42, 
          sizes = c(95, 48, 24, 12), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE48968-GPL13112sim123.json")

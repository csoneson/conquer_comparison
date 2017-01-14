suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE48968-GPL13112mock.rds", 
          subfile = "subsets/GSE48968-GPL13112mock_subsets.rds",
          resfilebase = "results/GSE48968-GPL13112mock",
          figfilebase = "figures/diffexpression/GSE48968-GPL13112mock", 
          groupid = "source_name_ch1", 
          keepgroups = "BMDC (1h LPS Stimulation)", 
          seed = 42, 
          sizes = c(48, 24, 12, 6), 
          nreps = c(1, 5, 5, 5))
write(toJSON(L), file = "config/GSE48968-GPL13112mock.json")

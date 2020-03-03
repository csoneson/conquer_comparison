suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE62270-GPL17021mock.rds", 
          subfile = "subsets/GSE62270-GPL17021mock_subsets.rds",
          resfilebase = "results/GSE62270-GPL17021mock",
          figfilebase = "figures/diffexpression/GSE62270-GPL17021mock", 
          groupid = "source_name_ch1", 
          keepgroups = c("Randomly extracted ex vivo isolated 5 day YFP positive cells"), 
          seed = 42, 
          sizes = c(200, 120, 48, 24, 12), 
          nreps = c(1, 5, 5, 5, 5))
write(toJSON(L), file = "config/GSE62270-GPL17021mock.json")

suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE41265.rds", 
          subfile = "subsets/GSE41265_subsets.rds",
          resfilebase = "results/GSE41265",
          figfilebase = "figures/GSE41265", 
          groupid = "source_name_ch1", 
          keepgroups = NULL, 
          seed = 42, 
          # sizes = c(50, 30, 24, 12), 
          # nreps = c(1, 3, 3, 3),
          sizes = c(4), 
          nreps = c(1))
write(toJSON(L), file = "config/GSE41265.json")


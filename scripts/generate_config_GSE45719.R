suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE45719.rds", 
          subfile = "subsets/GSE45719_subsets.rds",
          resfilebase = "results/GSE45719",
          figfilebase = "figures/GSE45719", 
          groupid = "source_name_ch1", 
          keepgroups = c("16-cell stage blastomere",
                         "Mid blastocyst cell (92-94h post-fertilization)"), 
          seed = 42, 
          # sizes = c(50, 30, 24, 12), 
          # nreps = c(1, 3, 3, 3),
          sizes = c(12), 
          nreps = c(1))
write(toJSON(L), file = "config/GSE45719.json")

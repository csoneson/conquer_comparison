suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE63818-GPL16791.rds", 
          subfile = "subsets/GSE63818-GPL16791_subsets.rds",
          resfilebase = "results/GSE63818-GPL16791",
          figfilebase = "figures/diffexpression/GSE63818-GPL16791", 
          groupid = c("source_name_ch1", "characteristics_ch1"), 
          keepgroups = c("Primordial Germ Cells.developmental stage: 7 week gestation",
                         "Somatic Cells.developmental stage: 7 week gestation"), 
          seed = 42, 
          sizes = c(26, 12, 6), 
          nreps = c(1, 5, 5))
write(toJSON(L), file = "config/GSE63818-GPL16791.json")

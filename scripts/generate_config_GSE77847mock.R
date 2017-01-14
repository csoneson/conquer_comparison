suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE77847mock.rds", 
          subfile = "subsets/GSE77847mock_subsets.rds",
          resfilebase = "results/GSE77847mock",
          figfilebase = "figures/diffexpression/GSE77847mock", 
          groupid = "characteristics_ch1", 
          keepgroups = "sample type: cKit+ Flt3ITD/ITD,Dnmt3afl/- MxCre  AML-1", 
          seed = 42, 
          sizes = c(24, 12, 6), 
          nreps = c(1, 5, 5))
write(toJSON(L), file = "config/GSE77847mock.json")

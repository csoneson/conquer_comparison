suppressPackageStartupMessages(library(rjson))

## Generate configuration file
L <- list(mae = "data/GSE77847.rds", 
          subfile = "subsets/GSE77847_subsets.rds",
          resfilebase = "results/GSE77847",
          figfilebase = "figures/diffexpression/GSE77847", 
          groupid = "characteristics_ch1", 
          keepgroups = c("sample type: cKit+ Flt3ITD/ITD,Dnmt3afl/- MxCre  AML-1",
                         "sample type: cKit+ Flt3ITD/ITD,Dnmt3afl/- MxCre  AML-2"), 
          seed = 42, 
          sizes = c(48, 36, 24, 12), 
          nreps = c(1, 3, 3, 3))
write(toJSON(L), file = "config/GSE77847.json")

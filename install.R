# CRAN packages
install.packages(
  c(
    "tidyverse",
    "devtools",
    "testthat",
    "rjson",
    "survey",
    "ggrepel",
    "Rtsne",
    "akima",
    "scatterplot3d"
  )
)

# install BioConductor packages
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(
  c(
    "Biostrings",
    "SummarizedExperiment",
    "MultiAssayExperiment",
    "iCOBRA",
    "monocle",
    "edgeR",
    "MAST",
    "limma",
    "genefilter",
    "scater",
    "IHW",
    "GEOquery"
  ))
  
# GitHub packages
devtools::install_github("vqv/ggbiplot")

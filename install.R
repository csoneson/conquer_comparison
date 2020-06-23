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
    "scatterplot3d",
    "statmod",
    "cowplot",
    "Matrix",
    "pheatmap",
    "reshape2",
    "cluster",
    "lazyeval",
    "samr",
    "UpSetR",
    "RColorBrewer",
    # Seurat dependencies
    "mixtools",
    "lars",
    "tsne",
    "fpc",
    "pbapply",
    "FNN",
    "caret",
    "RcppProgress",
    "tclust",
    "ranger"
  )
)

# install BioConductor packages
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(
  c(
    "Biobase",
    "Biostrings",
    "DESeq2",
    "SummarizedExperiment",
    "MultiAssayExperiment",
    "iCOBRA",
    "monocle",
    "edgeR",
    "MAST",
    "limma",
    "genefilter",
    "scran",
    "scater",
    "IHW",
    "GEOquery",
    "ggtree",
    "metagenomeSeq",
    "ROTS",
    "scde",
    "tximport",
    "impute"
  ))
  
# GitHub packages
devtools::install_github("vqv/ggbiplot")
devtools::install_github("nghiavtr/BPSC")

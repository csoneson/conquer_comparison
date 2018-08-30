## Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This repository contains all the necessary code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data, available in 

* C Soneson & MD Robinson: [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612). Nature Methods 15:255-261 (2018).

In this paper, we compare the performance of more than 30 approaches to differential gene expression analysis in the context of single-cell RNA-seq data. The main results can be further browsed in a [shiny app](http://imlspenticton.uzh.ch:3838/scrnaseq_de_evaluation). 

**Note:** The purpose of the `conquer_comparison` repository is to provide a public record of the exact code that was used for our publication ([Soneson & Robinson, Nature Methods 2018](https://www.nature.com/articles/nmeth.4612)). In particular, it is not intended to be a software package or a general pipeline for differential expression analysis of single-cell data. As a consequence, running the code requires the same software and package versions that were used for our analyses (all versions are indicated in the paper). As the analysis involved running a large number of methods on many data sets and over an extended period of time, we cannot guarantee that it will run successfully with new releases of the software, or that exactly the same results will be obtained with newer versions of the packages. While the repository will not be updated to ensure that it runs with every new version of the used packages, the issues can be used to post questions and/or solutions as they arise. 

The repository contains the following information:

* `config/` contains configuration files for all the data sets that we considered. The configuration files detail the cell populations that were compared, as well as the number of cells per group used in each comparison.
* `data/` contains some of the raw data that was used for the comparison. All data sets that were used can be downloaded as a bundle from [http://imlspenticton.uzh.ch/robinson_lab/conquer\_de\_comparison/](http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/)
* `export_results/` contains results for the final figures, in tabular format
* `scripts/` contains all R scripts used for the evaluation
* `shiny/` contains the code for a shiny app built to browse the results ([http://imlspenticton.uzh.ch:3838/scrnaseq\_de\_evaluation](http://imlspenticton.uzh.ch:3838/scrnaseq_de_evaluation))
* `unit_tests/` contains unit tests that were used to check the calculations
* `Makefile` is the master script, which outlines the entire evaluation and calls all scripts in the appropriate order
* `include_filterings.mk`, `include_datasets.mk`, `include_methods.mk` and `plot_methods.mk` are additional makefiles listing the filter settings, data set and differential expression methods used in the comparison 
 

## Running the comparison
Assuming that all prerequisites are available, the comparison can be run by simply typing 

```$ make```

from the top directory (note, however, that this will take a **significant** amount of time!). The Makefile reads the three files *include_filterings.mk*, *include_datasets.mk* and *include_methods.mk* and performs the evaluation using the data sets, methods and filterings defined in these. The file *plot_methods.mk* detail the methods included in the final summary plots. For the code to execute properly, an *.rds* file containing a *MultiAssayExperiment* object for each data set must be provided in the `data/` directory. Such files can be downloaded, e.g., from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) database. The files used for the evaluation are bundled together in an archive that can be downloaded from [here](http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/)

## Adding a differential expression method
To add a differential expression method to the evaluation, construct a script in the form of the provided `apply_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the method. Then add the name of the method to `include_methods.mk`. To make it show up in the summary plots, add it to `plot_methods.mk` and assign it a color in `scripts/plot_setup.R`.

## Adding a data set
To add a data set, put the *.rds* file containing the *MultiArrayExperiment* object in the `data/` folder and construct a script in the form of the provided `generate_config_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the data set. Then add the name of the dataset to the appropriate variables in `include_datasets.mk`. Also, add the data set to the `data/dataset_type.txt` file, indicating the type of values in each data set.

## A note on the data sets
Most data sets in the published evaluation are obtained from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) repository. The RPM values for the Usoskin dataset was downloaded from [http://linnarssonlab.org/drg/](http://linnarssonlab.org/drg/) on December 18, 2016. The 10X data set was downloaded from [https://support.10xgenomics.com/single-cell-gene-expression/datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets) on September 17, 2017.

## Cell cycle genes
The list of mouse cell cycle genes was obtained from [http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html](http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html) on March 9, 2017.

## Unit tests
To run all the unit tests, start `R`, load the `testthat` package and run 
``source("scripts/run_unit_tests.R")``. Alternatively, to run just the unit tests in a given file, do e.g. ``test_file("unit_tests/test_trueperformance.R", reporter = "summary")``.

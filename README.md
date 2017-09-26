## Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This repository contains all the necessary code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data, available in 

* C Soneson & MD Robinson: [Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data](http://biorxiv.org/content/early/2017/05/28/143289). 

In this paper, we compare the performance of 30 approaches to differential gene expression analysis in the context of single-cell RNA-seq data. 

The repository contains the following information:

* `config/` contains configuration files for all the data sets that we considered. The configuration files detail the cell populations that were compared, as well as the number of cells per group used in each comparison.
* `scripts/` contains all R scripts used for the evaluation
* `unit_tests/` contains unit tests that were used to check the calculations
* `Makefile` is the master script, which outlines the entire evaluation and calls all scripts in the appropriate order
* `include_filterings.mk`, `include_datasets.mk` and `include_methods.mk` are additional makefiles listing the filter settings, data set and differential expression methods used in the comparison 
 

## Running the comparison
Assuming that all prerequisites are available, the comparison can be run by simply typing 

```$ make```

from the top directory (note, however, that this will take a **significant** amount of time!). The Makefile reads the three files *include_filterings.mk*, *include_datasets.mk* and *include_methods.mk* and performs the evaluation using the data sets, methods and filterings defined in these. For the code to execute properly, an *.rds* file containing a *MultiAssayExperiment* object for each data set must be provided in the `data/` directory. Such files can be downloaded, e.g., from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) database. 

## Adding a differential expression method
To add a differential expression method to the evaluation, construct a script in the form of the provided `apply_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the method. Then add the name of the method to `include_methods.mk`. To make it show up in the summary plots, assign it a color in `scripts/plot_setup.R`.

## Adding a data set
To add a data set, put the *.rds* file containing the *MultiArrayExperiment* object in the `data/` folder and construct a script in the form of the provided `generate_config_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the data set. Then add the name of the dataset to the appropriate variables in `include_datasets.mk`. Also, add the data set to the `data/dataset_type.txt` file, indicating the type of values in each data set.

## A note on the data sets
Most data sets in the published evaluation are obtained from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) repository. The RPM values for the Usoskin dataset was downloaded from [http://linnarssonlab.org/drg/](http://linnarssonlab.org/drg/) on December 18, 2016. 

## Cell cycle genes
The list of mouse cell cycle genes was obtained from [http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html](http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html) on March 9, 2017.

## Unit tests
To run all the unit tests, start `R`, load the `testthat` package and run 
``source("scripts/run_unit_tests.R")``. Alternatively, to run just the unit tests in a given file, do e.g. ``test_file("unit_tests/test_trueperformance.R", reporter = "summary")``.

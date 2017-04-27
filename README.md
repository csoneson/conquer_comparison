This repository contains the code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data. 

## Running the comparison
To run the comparison and generate plots, just type 

```$ make```

from the top directory. The Makefile reads the three files *include_filterings.mk*, *include_datasets.mk* and *include_methods.mk* and performs the evaluation using the data sets, methods and filterings defined in these. For the code to execute, an *.rds* file containing a *MultiAssayExperiment* object for each data set must be provided in the `data/` directory. Such files can be downloaded, e.g., from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) database. 

## Adding a method
To add a method to the evaluation, construct a script in the form of the provided `apply_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the method. Then add the name of the method to `include_methods.mk`.

## Adding a dataset
To add a dataset, put the *.rds* file containing the *MultiArrayExperiment* object in the `data/` folder and construct a script in the form of the provided `generate_config_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the dataset. Then add the name of the dataset to `include_datasets.mk`. 

## A note on the datasets
Most datasets in the evaluation are obtained from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) repository. The RPM values for the Usoskin dataset was downloaded from [http://linnarssonlab.org/drg/](http://linnarssonlab.org/drg/) on December 18, 2016. 

## Unit tests
To run all the unit tests, just do (in R)
``source("scripts/run_unit_tests.R")``. Alternatively, to run just the unit tests in a given file, do e.g. ``test_file("unit_tests/test_trueperformance.R", reporter = "summary")``.
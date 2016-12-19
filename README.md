## Running the comparison
To run the comparison and generate plots with the methods and datasets defined in `include_methods.mk`, just type 

```$ make```

## Adding a method
To add a method, construct a script in the form of the `apply_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the method. Then add the name of the method to `include_methods.mk`

## Adding a dataset
To add a dataset, put the rds file containing the MultiArrayExperiment in the `data/` folder and construct a script in the form of the `generate_config_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the dataset. Then add the name of the dataset to `include_methods.mk`. Note that this implies that (for simplicity), when a dataset is added, the plots for all the datasets in `include_methods.mk` will be regenerated. 

## A note on the datasets
Most datasets are obtained from the [*conquer*](http://imlspenticton.uzh.ch:3838/conquer/) repository. The RPM values for the Usoskin dataset was downloaded from [http://linnarssonlab.org/drg/](http://linnarssonlab.org/drg/) on December 18, 2016. 
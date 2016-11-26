## Define the version of R and the path to the library
R := R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define the data sets to run the analysis for
#DS := GSE45719 

## Define the methods to apply to the data set
#MT := edgeRLRT SAMseq Wilcoxon edgeRZILRT edgeRQLF NODES monocle monoclecounts BPSC DESeq2 edgeRLRTdeconv edgeRLRTrobust MASTcounts MASTcountsDetRate MASTtpm SCDE 
#comma := ,
#empty :=
#space := $(empty) $(empty)
#MTc := $(subst $(space),$(comma),$(MT))
include include_methods.mk

.PHONY: all

## Define the default rule
#all: $(addsuffix .rds, $(addprefix results/, $(foreach X,$(DS),$(foreach Y,$(MT),$X_$Y))))
all: $(addsuffix .pdf, $(addprefix figures/comparison/, $(foreach X,$(DS),$X))) 

## Make sure no intermediate files are deleted
.SECONDARY:

## Generate configuration files
config/%.json: scripts/generate_config_%.R
	$R scripts/generate_config_$*.R Rout/generate_config_$*.Rout

## Extract sample subsets
subsets/%_subsets.rds: data/%.rds config/%.json scripts/generate_subsets.R
	$R "--args config_file='config/$*.json'" scripts/generate_subsets.R Rout/generate_subsets_$*.Rout

## Define rules for differential expression
define dgerule
results/$(1)_$(2).rds: scripts/apply_$(2).R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds
	$R "--args config_file='config/$(1).json' demethod='$(2)'" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$(2).Rout
endef
$(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerule,$(i),$(j)))))

## Plots for comparison
figures/comparison/%.pdf: $(addsuffix .rds, $(addprefix results/%_, $(foreach Y,$(MT),$Y))) scripts/plot_comparison.R scripts/plot_functions.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$*' config_file='config/$*.json'" scripts/plot_comparison.R Rout/plot_comparison_$*.Rout

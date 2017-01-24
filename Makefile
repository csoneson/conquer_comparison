## Define the version of R and the path to the library
R := R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define the active datasets and methods
include include_methods.mk

.PHONY: all

## Define the default rule
all: $(addsuffix .pdf, $(addprefix figures/comparison/, $(foreach X,$(DS),$X))) \
$(addsuffix .pdf, $(addprefix figures/comparison/, $(foreach Y,$(FILT),$(foreach X,$(DS),$X_$Y)))) \
$(addsuffix _orig_vs_mock.pdf, $(addprefix figures/orig_vs_mock/, $(foreach X,$(DSb),$X))) \
$(addsuffix _orig_vs_mock.pdf, $(addprefix figures/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(DSb),$X_$Y)))) \
figures/summary_crossds/summary_heatmaps.pdf \
$(addsuffix .pdf, $(addprefix figures/dataset_characteristics/, $(foreach X,$(DS),$X))) \
$(addsuffix .pdf, $(addprefix figures/dataset_characteristics/, $(foreach Y,$(FILT),$(foreach X,$(DS),$X_$Y))))

## Make sure no intermediate files are deleted
.SECONDARY:

## Generate configuration files
config/%.json: scripts/generate_config_%.R
	$R scripts/generate_config_$*.R Rout/generate_config_$*.Rout

## Extract sample subsets
subsets/%_subsets.rds: data/%.rds config/%.json scripts/generate_subsets.R
	$R "--args config_file='config/$*.json'" scripts/generate_subsets.R Rout/generate_subsets_$*.Rout

## Generate Usoskin data set
data/UsoskinGSE59739.rds: scripts/generate_Usoskin_mae.R data/Usoskin_External_resources_Table_1.txt
	$R scripts/generate_Usoskin_mae.R Rout/generate_Usoskin_mae.Rout

## Define rules for differential expression
## Without filtering
define dgerule
results/$(1)_$(2).rds: scripts/apply_$(2).R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args config_file='config/$(1).json' demethod='$(2)' filt=''" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$(2).Rout
endef
$(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerule,$(i),$(j)))))

## With filtering
define dgerulefilt
results/$(1)_$(2)_$(3).rds: scripts/apply_$(2).R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args config_file='config/$(1).json' demethod='$(2)' filt='$(3)'" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$(2)_$(3).Rout
endef
$(foreach k, $(FILT), $(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerulefilt,$(i),$(j),$(k))))))

## Plots for comparison
figures/comparison/%.pdf: $(addsuffix .rds, $(addprefix results/%_, $(foreach Y,$(MT),$Y))) \
scripts/plot_comparison.R scripts/plot_functions.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$*' config_file='config/$*.json' filt=''" scripts/plot_comparison.R Rout/plot_comparison_$*.Rout

define plotrule
figures/comparison/$(1)_$(2).pdf: $(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) \
scripts/plot_comparison.R scripts/plot_functions.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$(1)' config_file='config/$(1).json' filt='$(2)'" scripts/plot_comparison.R Rout/plot_comparison_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(DS),$(eval $(call plotrule,$(i),$(k)))))

## Plots for characterization of data set
figures/dataset_characteristics/%.pdf: include_methods.mk scripts/plot_characterize_dataset.R scripts/prepare_mae.R \
subsets/%_subsets.rds data/%.rds
	$R "--args dataset='$*' config_file='config/$*.json' filt=''" scripts/plot_characterize_dataset.R Rout/plot_characterize_dataset_$*.Rout

define plotrule_characterization
figures/dataset_characteristics/$(1)_$(2).pdf: include_methods.mk scripts/plot_characterize_dataset.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args dataset='$(1)' config_file='config/$*.json' filt='$(2)'" scripts/plot_characterize_dataset.R Rout/plot_characterize_dataset_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(DS),$(eval $(call plotrule_characterization,$(i),$(k)))))

## Plots for comparison, orig vs mock
figures/orig_vs_mock/%_orig_vs_mock.pdf: $(addsuffix .rds, $(addprefix results/%_, $(foreach Y,$(MT),$Y))) \
$(addsuffix .rds, $(addprefix results/%mock_, $(foreach Y,$(MT),$Y))) \
scripts/plot_orig_vs_mock.R scripts/plot_functions.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$*' filt=''" scripts/plot_orig_vs_mock.R Rout/plot_orig_vs_mock_$*.Rout

define plotrule_origvsmock
figures/orig_vs_mock/$(1)_$(2)_orig_vs_mock.pdf: $(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) \
$(addsuffix _$(2).rds, $(addprefix results/$(1)mock_, $(foreach Y,$(MT),$Y))) \
scripts/plot_orig_vs_mock.R scripts/plot_functions.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$(1)' filt='$(2)'" scripts/plot_orig_vs_mock.R Rout/plot_comparison_orig_vs_mock_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(DSb),$(eval $(call plotrule_origvsmock,$(i),$(k)))))

## Plots for comparison of characteristics of significant genes, across mock data sets
figures/summary_crossds/summary_heatmaps.pdf: $(addsuffix .pdf, $(addprefix figures/comparison/, $(foreach Y,$(Dss),$Y))) \
scripts/plot_summarize_datasets.R include_methods.mk
	$R "--args datasets='${Dssc}' filt=''" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets.Rout




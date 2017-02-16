## Define the version of R and the path to the library
R := R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define the active datasets and methods
include include_methods.mk

## Plot types
PLOTTYPE := ks timing truefpr results_characterization consistency results_relativetruth results_relativetruth_all
SUMMARYTYPE := truefpr pca timing fracNA

.PHONY: all

## Define the default rule
all: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach k,$(PLOTTYPE),$(foreach X,$(DS),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach k,$(PLOTTYPE),$(foreach Y,$(FILT),$(foreach X,$(DS),$k/$X_$Y_$k))))) \
$(addsuffix .pdf, $(addprefix figures/dataset_characteristics/, $(foreach X,$(DS),$X))) \
$(addsuffix .pdf, $(addprefix figures/dataset_characteristics/, $(foreach Y,$(FILT),$(foreach X,$(DS),$X_$Y)))) \
$(addsuffix .rds, $(addprefix figures/summary_crossds/summary_, $(foreach K,$(SUMMARYTYPE),$(K)))) \
$(addsuffix .rds, $(addprefix figures/summary_crossds/summary_, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)_$(Y))))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix figures/orig_vs_mock/, $(foreach X,$(Dsb),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix figures/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(Dsb),$X_$Y)))) \
figures/summary_crossds/summary_orig_vs_mock.rds \
$(addsuffix _orig_vs_mock.rds, $(addprefix figures/summary_crossds/summary_, $(foreach Y,$(FILT),$Y)))

## Make sure no intermediate files are deleted
.SECONDARY:

## -------------------------- Generate configuration files ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
config/%.json: scripts/generate_config_%.R
	$R scripts/generate_config_$*.R Rout/generate_config_$*.Rout

## --------------------------- Extract sample subsets --------------------------------- ##
## ------------------------------------------------------------------------------------ ##
subsets/%_subsets.rds: data/%.rds config/%.json scripts/generate_subsets.R
	$R "--args config_file='config/$*.json'" scripts/generate_subsets.R Rout/generate_subsets_$*.Rout

## -------------------------- Generate Usoskin data set ------------------------------- ##
## ------------------------------------------------------------------------------------ ##
data/UsoskinGSE59739.rds: scripts/generate_Usoskin_mae.R data/Usoskin_External_resources_Table_1.txt
	$R scripts/generate_Usoskin_mae.R Rout/generate_Usoskin_mae.Rout

## ------------------ Define rules for differential expression ------------------------ ##
## ------------------------------------------------------------------------------------ ##
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

## ------------------ Prepare COBRAData object for evaluation ------------------------- ##
## ------------------------------------------------------------------------------------ ##
figures/cobra_data/%_cobra.rds: include_methods.mk scripts/prepare_cobra_for_evaluation.R \
$(addsuffix .rds, $(addprefix results/%_, $(foreach Y,$(MT),$Y))) scripts/prepare_mae.R
	$R "--args demethods='${MTc}' dataset='$*' config_file='config/$*.json' filt=''" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$*.Rout

define cobrarule_filt
figures/cobra_data/$(1)_$(2)_cobra.rds: include_methods.mk scripts/prepare_cobra_for_evaluation.R \
$(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) scripts/prepare_mae.R
	$R "--args demethods='${MTc}' dataset='$(1)' config_file='config/$(1).json' filt='$(2)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call cobrarule_filt,$(X),$(k)))))

## --------------------------- Plots for evaluation ----------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define plotrule
figures/$(2)/$(1)_$(2)_summary_data.rds: scripts/plot_evaluation.R scripts/plot_$(2).R scripts/plot_setup.R figures/cobra_data/$(1)_cobra.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' plottype='$(2)'" scripts/plot_evaluation.R Rout/plot_evaluation_$(1)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE),$(eval $(call plotrule,$(X),$(Y)))))

define plotrule_filt
figures/$(2)/$(1)_$(3)_$(2)_summary_data.rds: scripts/plot_evaluation.R scripts/plot_$(2).R scripts/plot_setup.R figures/cobra_data/$(1)_$(3)_cobra.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(3)' plottype='$(2)'" scripts/plot_evaluation.R Rout/plot_evaluation_$(1)_$(3)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE),$(eval $(call plotrule_filt,$(X),$(Y),$(k))))))

## -------------------- Plots for characterization of data set ------------------------ ##
## ------------------------------------------------------------------------------------ ##
define plotrule_characterization
figures/dataset_characteristics/$(1).pdf: include_methods.mk scripts/plot_characterize_dataset.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt=''" scripts/plot_characterize_dataset.R Rout/plot_characterize_dataset_$(1).Rout
endef
$(foreach i,$(DS),$(eval $(call plotrule_characterization,$(i))))

define plotrule_characterization_filt
figures/dataset_characteristics/$(1)_$(2).pdf: include_methods.mk scripts/plot_characterize_dataset.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(2)'" scripts/plot_characterize_dataset.R Rout/plot_characterize_dataset_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(DS),$(eval $(call plotrule_characterization_filt,$(i),$(k)))))

## -------------------- Plots for evaluation, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
figures/orig_vs_mock/%_orig_vs_mock_summary_data.rds: $(addsuffix .rds, $(addprefix results/%_, $(foreach Y,$(MT),$Y))) \
$(addsuffix .rds, $(addprefix results/%mock_, $(foreach Y,$(MT),$Y))) \
scripts/plot_orig_vs_mock.R 
	$R "--args demethods='${MTc}' dataset='$*' filt=''" scripts/plot_orig_vs_mock.R Rout/plot_orig_vs_mock_$*.Rout

define plotrule_origvsmock
figures/orig_vs_mock/$(1)_$(2)_orig_vs_mock_summary_data.rds: $(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) \
$(addsuffix _$(2).rds, $(addprefix results/$(1)mock_, $(foreach Y,$(MT),$Y))) \
scripts/plot_orig_vs_mock.R 
	$R "--args demethods='${MTc}' dataset='$(1)' filt='$(2)'" scripts/plot_orig_vs_mock.R Rout/plot_orig_vs_mock_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(Dsb),$(eval $(call plotrule_origvsmock,$(i),$(k)))))

## ---------------------- Summary plots, across mock data sets ------------------------ ##
## ------------------------------------------------------------------------------------ ##
figures/summary_crossds/summary_truefpr.rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),truefpr/$Y_truefpr))) \
scripts/plot_summarize_datasets.R scripts/summarize_truefpr.R
	$R "--args datasets='${Dssc}' filt='' summarytype='truefpr'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_truefpr.Rout

figures/summary_crossds/summary_pca.rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),results_characterization/$Y_results_characterization))) \
scripts/plot_summarize_datasets.R scripts/summarize_pca.R
	$R "--args datasets='${Dssc}' filt='' summarytype='pca'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_pca.Rout

figures/summary_crossds/summary_timing.rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),timing/$Y_timing))) \
scripts/plot_summarize_datasets.R scripts/summarize_timing.R
	$R "--args datasets='${Dssc}' filt='' summarytype='timing'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_timing.Rout

figures/summary_crossds/summary_fracNA.rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),cobra_data/$Y))) \
scripts/plot_summarize_datasets.R scripts/summarize_fracNA.R
	$R "--args datasets='${Dssc}' filt='' summarytype='fracNA'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_fracNA.Rout

define summaryrule_truefpr
figures/summary_crossds/summary_truefpr_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),truefpr/$Y_$(1)_truefpr))) \
scripts/plot_summarize_datasets.R scripts/summarize_truefpr.R
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='truefpr'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_truefpr_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr,$(k))))

define summaryrule_pca
figures/summary_crossds/summary_pca_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),results_characterization/$Y_$(1)_results_characterization))) \
scripts/plot_summarize_datasets.R scripts/summarize_pca.R
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='pca' summarytype='pca'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_pca_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_pca,$(k))))

define summaryrule_timing
figures/summary_crossds/summary_timing_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),timing/$Y_$(1)_timing))) \
scripts/plot_summarize_datasets.R scripts/summarize_timing.R
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='timing' summarytype='timing'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_timing_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_timing,$(k))))

define summaryrule_fracna
figures/summary_crossds/summary_fracNA_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix figures/, $(foreach Y,$(Dss),cobra_data/$Y_$(1)))) \
scripts/plot_summarize_datasets.R scripts/summarize_fracNA.R
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='fracNA' summarytype='fracNA'" scripts/plot_summarize_datasets.R Rout/plot_summarize_datasets_fracNA_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_fracna,$(k))))

## --------------------------- Summary plots, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
figures/summary_crossds/summary_orig_vs_mock.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix figures/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/plot_summarize_orig_vs_mock.R
	$R "--args datasets='${Dsbc}' filt=''" scripts/plot_summarize_orig_vs_mock.R Rout/plot_summarize_orig_vs_mock.Rout

define plotrule_summary_origvsmock
figures/summary_crossds/summary_$(1)_orig_vs_mock.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix figures/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/plot_summarize_orig_vs_mock.R
	$R "--args datasets='${Dsbc}' filt='$(1)'" scripts/plot_summarize_orig_vs_mock.R Rout/plot_summarize_$(1)_orig_vs_mock.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmock,$(k))))



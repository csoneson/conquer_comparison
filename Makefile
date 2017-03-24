## Define the version of R and the path to the library
R := R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define paths
cobradir := output/cobra_data
concordancedir := output/concordances
relperfdir := output/relative_performance
realperfdir := output/performance_realtruth
figdir := figures
singledsfigdir := figures/single_dataset
multidsfigdir := figures/multi_dataset
dschardir := figures/dataset_characteristics

## Define the active datasets and methods
include include_methods.mk
include include_datasets.mk
include include_filterings.mk

## Plot types
PLOTTYPE := timing results_characterization results_relativetruth consistency results_relativetruth_all truefpr
PLOTTYPE1 := timing results_characterization results_relativetruth truefpr
PLOTTYPE2 := consistency
PLOTTYPE3 := results_relativetruth_all
PLOTTYPE4 := performance_realtruth
SUMMARYTYPE := truefpr pca timing fracNA crossmethod_consistency relfprtpr
SUMMARYTYPE2 := filtering
SUMMARYTYPE3 := trueperformance

.PHONY: all

## Define the default rule
all: plotds plotind plotorigmock \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock.rds \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_sim.rds \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_bulk.rds \
$(addsuffix .rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix .rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix .rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _bulk.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix _bulk.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _bulk.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix .rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _bulk.rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE3),$(K)/summary_$(K)))) \
$(addsuffix _sim.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE3),$(K)/summary_$(K)_$(Y)))))

## Plot original vs mock comparison
plotorigmock: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(Dsb),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(Dsb),$X_$Y)))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(Dsbsim),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(Dsbsim),$X_$Y)))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(DSbulkb),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(DSbulkb),$X_$Y))))

## Plot individual data set results
plotind: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach X,$(DS),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach Y,$(FILT),$(foreach X,$(DS),$k/$X_$Y_$k))))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE4),$(foreach X,$(DSsim),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE4),$(foreach Y,$(FILT),$(foreach X,$(DSsim),$k/$X_$Y_$k))))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach X,$(DSbulk),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach Y,$(FILT),$(foreach X,$(DSbulk),$k/$X_$Y_$k)))))

## Plot data set characteristics
plotds: $(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach X,$(DS),$X))) \
$(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(FILT),$(foreach X,$(DS),$X_$Y)))) \
$(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach X,$(DSbulk),$X))) \
$(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(FILT),$(foreach X,$(DSbulk),$X_$Y))))

## Prepare results for plotting, step II
plotprepare: $(addsuffix _concordances.rds, $(addprefix $(concordancedir)/, $(foreach k,$(DS),$k))) \
$(addsuffix _concordances.rds, $(addprefix $(concordancedir)/, $(foreach j,$(FILT),$(foreach k,$(DS),$k_$j)))) \
$(addsuffix _concordances.rds, $(addprefix $(concordancedir)/, $(foreach k,$(DSbulk),$k))) \
$(addsuffix _concordances.rds, $(addprefix $(concordancedir)/, $(foreach j,$(FILT),$(foreach k,$(DSbulk),$k_$j)))) \
$(addsuffix _relative_performance.rds, $(addprefix $(relperfdir)/, $(foreach k,$(DS),$k))) \
$(addsuffix _relative_performance.rds, $(addprefix $(relperfdir)/, $(foreach j,$(FILT),$(foreach k,$(DS),$k_$j)))) \
$(addsuffix _relative_performance.rds, $(addprefix $(relperfdir)/, $(foreach k,$(DSbulk),$k))) \
$(addsuffix _relative_performance.rds, $(addprefix $(relperfdir)/, $(foreach j,$(FILT),$(foreach k,$(DSbulk),$k_$j)))) \
$(addsuffix _performance.rds, $(addprefix $(realperfdir)/, $(foreach k,$(DSsim),$k))) \
$(addsuffix _performance.rds, $(addprefix $(realperfdir)/, $(foreach j,$(FILT),$(foreach k,$(DSsim),$k_$j))))

## Prepare results for plotting, step I
cobra: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach k,$(DS),$k))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach j,$(FILT),$(foreach k,$(DS),$k_$j)))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach k,$(DSbulk),$k))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach j,$(FILT),$(foreach k,$(DSbulk),$k_$j))))

## Run differential expression
diffexp: $(addsuffix .rds, $(addprefix results/, $(foreach k,$(MT),$(foreach X,$(DS),$X_$k)))) \
$(addsuffix .rds, $(addprefix results/, $(foreach Y,$(FILT),$(foreach k,$(MT),$(foreach X,$(DS),$X_$k_$Y))))) \
$(addsuffix .rds, $(addprefix results/, $(foreach k,$(MTbulk),$(foreach X,$(DSbulk),$X_$k)))) \
$(addsuffix .rds, $(addprefix results/, $(foreach Y,$(FILT),$(foreach k,$(MTbulk),$(foreach X,$(DSbulk),$X_$k_$Y)))))

## Simulate data
sim: $(addsuffix .rds, $(addprefix data/, $(foreach X,$(DSforsim),$Xsim123))) \
$(addsuffix .rds, $(addprefix data/, $(foreach X,$(DSforsim),$Xsim123mock)))

## List all rules
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

## Make sure no intermediate files are deleted
.SECONDARY:

## -------------------------- Generate configuration files ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
define configrule
config/$(1).json: scripts/generate_config_$(1).R
	$R scripts/generate_config_$(1).R Rout/generate_config_$(1).Rout
endef
$(foreach j,$(DS),$(eval $(call configrule,$(j))))
$(foreach j,$(DSbulk),$(eval $(call configrule,$(j))))

## --------------------------------- Simulate data ------------------------------------ ##
## ------------------------------------------------------------------------------------ ##
define simrule
data/$(1)sim$(2).rds: scripts/simulate_data.R data/$(1).rds config/$(1).json scripts/powsim_modified_functions.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' pDE=0.1 seed=$(2)" scripts/simulate_data.R Rout/simulate_data_$(1)_$(2).Rout
endef
$(foreach j,$(DSforsim),$(eval $(call simrule,$(j),123)))

define simrulemock
data/$(1)sim$(2)mock.rds: data/$(1)sim$(2).rds
	scp data/$(1)sim$(2).rds data/$(1)sim$(2)mock.rds
endef
$(foreach j,$(DSforsim),$(eval $(call simrulemock,$(j),123)))

## --------------------------- Extract sample subsets --------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define subsetrule
subsets/$(1)_subsets.rds: data/$(1).rds config/$(1).json scripts/generate_subsets.R
	$R "--args config_file='config/$(1).json'" scripts/generate_subsets.R Rout/generate_subsets_$(1).Rout
endef
$(foreach j,$(DS),$(eval $(call subsetrule,$(j))))
$(foreach j,$(DSbulk),$(eval $(call subsetrule,$(j))))

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
$(foreach j,$(MTbulk), $(foreach i,$(DSbulk),$(eval $(call dgerule,$(i),$(j)))))

## With filtering
define dgerulefilt
results/$(1)_$(2)_$(3).rds: scripts/apply_$(2).R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args config_file='config/$(1).json' demethod='$(2)' filt='$(3)'" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$(2)_$(3).Rout
endef
$(foreach k, $(FILT), $(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerulefilt,$(i),$(j),$(k))))))
$(foreach k, $(FILT), $(foreach j,$(MTbulk), $(foreach i,$(DSbulk),$(eval $(call dgerulefilt,$(i),$(j),$(k))))))

## ------------------ Prepare COBRAData object for evaluation ------------------------- ##
## ------------------------------------------------------------------------------------ ##
define cobrarule
$(cobradir)/$(1)_cobra.rds: scripts/prepare_cobra_for_evaluation.R scripts/calculate_gene_characteristics.R \
$(addsuffix .rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) scripts/prepare_mae.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$(1)' config_file='config/$(1).json' filt='' resdir='results' outdir='$(cobradir)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1).Rout
endef
$(foreach X,$(DS),$(eval $(call cobrarule,$(X))))

define cobrarulebulk
$(cobradir)/$(1)_cobra.rds: scripts/prepare_cobra_for_evaluation.R scripts/calculate_gene_characteristics.R \
$(addsuffix .rds, $(addprefix results/$(1)_, $(foreach Y,$(MTbulk),$Y))) scripts/prepare_mae.R include_methods.mk
	$R "--args demethods='${MTcbulk}' dataset='$(1)' config_file='config/$(1).json' filt='' resdir='results' outdir='$(cobradir)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1).Rout
endef
$(foreach X,$(DSbulk),$(eval $(call cobrarulebulk,$(X))))

define cobrarule_filt
$(cobradir)/$(1)_$(2)_cobra.rds: scripts/prepare_cobra_for_evaluation.R scripts/calculate_gene_characteristics.R \
$(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MT),$Y))) scripts/prepare_mae.R include_methods.mk
	$R "--args demethods='${MTc}' dataset='$(1)' config_file='config/$(1).json' filt='$(2)' resdir='results' outdir='$(cobradir)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call cobrarule_filt,$(X),$(k)))))

define cobrarulebulk_filt
$(cobradir)/$(1)_$(2)_cobra.rds: scripts/prepare_cobra_for_evaluation.R scripts/calculate_gene_characteristics.R \
$(addsuffix _$(2).rds, $(addprefix results/$(1)_, $(foreach Y,$(MTbulk),$Y))) scripts/prepare_mae.R include_methods.mk
	$R "--args demethods='${MTcbulk}' dataset='$(1)' config_file='config/$(1).json' filt='$(2)' resdir='results' outdir='$(cobradir)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call cobrarulebulk_filt,$(X),$(k)))))

## ----------------------------- Calculate concordances ------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define concrule
$(concordancedir)/$(1)_concordances.rds: scripts/calculate_concordances.R $(cobradir)/$(1)_cobra.rds
	$R "--args dataset='$(1)' filt='' cobradir='$(cobradir)' outdir='$(concordancedir)'" scripts/calculate_concordances.R Rout/calculate_concordances_$(1).Rout
endef
$(foreach X,$(DS),$(eval $(call concrule,$(X))))
$(foreach X,$(DSbulk),$(eval $(call concrule,$(X))))

define concrule_filt
$(concordancedir)/$(1)_$(2)_concordances.rds: scripts/calculate_concordances.R $(cobradir)/$(1)_$(2)_cobra.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(concordancedir)'" scripts/calculate_concordances.R Rout/calculate_concordances_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call concrule_filt,$(X),$(k)))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call concrule_filt,$(X),$(k)))))

## ------------------------ Calculate relative performances --------------------------- ##
## ------------------------------------------------------------------------------------ ##
define relperfrule
$(relperfdir)/$(1)_relative_performance.rds: scripts/calculate_relative_performance_all_truths.R $(cobradir)/$(1)_cobra.rds
	$R "--args dataset='$(1)' filt='' cobradir='$(cobradir)' outdir='$(relperfdir)'" scripts/calculate_relative_performance_all_truths.R Rout/calculate_relative_performance_all_truths_$(1).Rout
endef
$(foreach X,$(DS),$(eval $(call relperfrule,$(X))))
$(foreach X,$(DSbulk),$(eval $(call relperfrule,$(X))))

define relperfrule_filt
$(relperfdir)/$(1)_$(2)_relative_performance.rds: scripts/calculate_relative_performance_all_truths.R $(cobradir)/$(1)_$(2)_cobra.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(relperfdir)'" scripts/calculate_relative_performance_all_truths.R Rout/calculate_relative_performance_all_truths_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call relperfrule_filt,$(X),$(k)))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call relperfrule_filt,$(X),$(k)))))

## -------------------------- Calculate true performances ----------------------------- ##
## ------------------------------------------------------------------------------------ ##
define trueperfrule
$(realperfdir)/$(1)_performance.rds: scripts/calculate_performance_realtruth.R $(cobradir)/$(1)_cobra.rds data/$(1)_truth.rds
	$R "--args dataset='$(1)' filt='' cobradir='$(cobradir)' outdir='$(realperfdir)'" scripts/calculate_performance_realtruth.R Rout/calculate_performance_realtruth_$(1).Rout
endef
$(foreach X,$(DSsim),$(eval $(call trueperfrule,$(X))))

define trueperfrule_filt
$(realperfdir)/$(1)_$(2)_performance.rds: scripts/calculate_performance_realtruth.R $(cobradir)/$(1)_$(2)_cobra.rds data/$(1)_truth.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(realperfdir)'" scripts/calculate_performance_realtruth.R Rout/calculate_performance_realtruth_$(1)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DSsim),$(eval $(call trueperfrule_filt,$(X),$(k)))))

## --------------------------- Plots for evaluation ----------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define plotrule
$(singledsfigdir)/$(2)/$(1)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)_cobra.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule,$(X),$(Y)))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule,$(X),$(Y)))))

define plotrule2
$(singledsfigdir)/$(2)/$(1)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)_cobra.rds \
$(concordancedir)/$(1)_concordances.rds scripts/help_function_crossmethod_concordance.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y)))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y)))))

define plotrule3
$(singledsfigdir)/$(2)/$(1)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(relperfdir)/$(1)_relative_performance.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y)))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y)))))

define plotrule4
$(singledsfigdir)/$(2)/$(1)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(realperfdir)/$(1)_performance.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(2).Rout
endef
$(foreach X,$(DSsim),$(foreach Y,$(PLOTTYPE4),$(eval $(call plotrule4,$(X),$(Y)))))

define plotrule_filt
$(singledsfigdir)/$(2)/$(1)_$(3)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)_$(3)_cobra.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(3)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule_filt,$(X),$(Y),$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule_filt,$(X),$(Y),$(k))))))

define plotrule2_filt
$(singledsfigdir)/$(2)/$(1)_$(3)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)_$(3)_cobra.rds \
$(concordancedir)/$(1)_$(3)_concordances.rds scripts/help_function_crossmethod_concordance.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(3)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2_filt,$(X),$(Y),$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2_filt,$(X),$(Y),$(k))))))

define plotrule3_filt
$(singledsfigdir)/$(2)/$(1)_$(3)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(relperfdir)/$(1)_$(3)_relative_performance.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(3)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3_filt,$(X),$(Y),$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3_filt,$(X),$(Y),$(k))))))

define plotrule4_filt
$(singledsfigdir)/$(2)/$(1)_$(3)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(realperfdir)/$(1)_$(3)_performance.rds
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)_$(3)_$(2).Rout
endef
$(foreach k,$(FILT),$(foreach X,$(DSsim),$(foreach Y,$(PLOTTYPE4),$(eval $(call plotrule4_filt,$(X),$(Y),$(k))))))

## -------------------- Plots for characterization of data set ------------------------ ##
## ------------------------------------------------------------------------------------ ##
define plotrule_characterization
$(dschardir)/$(1)_dataset_characteristics_summary_data.rds: scripts/run_plot_dataset_characterization.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds scripts/calculate_gene_characteristics.R  scripts/calculate_cell_characteristics.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='' cell_cycle_file='data/cell_cycle_geneids.rds' figdir='$(dschardir)'" scripts/run_plot_dataset_characterization.R Rout/run_plot_dataset_characterization_$(1).Rout
endef
$(foreach i,$(DS),$(eval $(call plotrule_characterization,$(i))))
$(foreach i,$(DSbulk),$(eval $(call plotrule_characterization,$(i))))

define plotrule_characterization_filt
$(dschardir)/$(1)_$(2)_dataset_characteristics_summary_data.rds: scripts/run_plot_dataset_characterization.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds scripts/calculate_gene_characteristics.R  scripts/calculate_cell_characteristics.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(2)' cell_cycle_file='data/cell_cycle_geneids.rds' figdir='$(dschardir)'" scripts/run_plot_dataset_characterization.R Rout/run_plot_dataset_characterization_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(DS),$(eval $(call plotrule_characterization_filt,$(i),$(k)))))
$(foreach k,$(FILT), $(foreach i,$(DSbulk),$(eval $(call plotrule_characterization_filt,$(i),$(k)))))

## -------------------- Plots for evaluation, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
define origvsmockrule
$(figdir)/orig_vs_mock/$(1)_orig_vs_mock_summary_data.rds: $(concordancedir)/$(1)_concordances.rds $(concordancedir)/$(1)mock_concordances.rds \
scripts/plot_setup.R scripts/run_plot_single_dataset_origvsmock.R scripts/plot_single_dataset_origvsmock.R
	$R "--args dataset='$(1)' filt='' concordancedir='$(concordancedir)' figdir='$(figdir)/orig_vs_mock'" scripts/run_plot_single_dataset_origvsmock.R Rout/run_plot_single_dataset_origvsmock_$(1).Rout
endef
$(foreach i,$(Dsb),$(eval $(call origvsmockrule,$(i))))
$(foreach i,$(DSbulkb),$(eval $(call origvsmockrule,$(i))))
$(foreach i,$(Dsbsim),$(eval $(call origvsmockrule,$(i))))

define origvsmockrule_filt
$(figdir)/orig_vs_mock/$(1)_$(2)_orig_vs_mock_summary_data.rds: $(concordancedir)/$(1)_$(2)_concordances.rds \
$(concordancedir)/$(1)mock_$(2)_concordances.rds scripts/plot_setup.R  \
scripts/run_plot_single_dataset_origvsmock.R scripts/plot_single_dataset_origvsmock.R
	$R "--args dataset='$(1)' filt='$(2)' concordancedir='$(concordancedir)' figdir='$(figdir)/orig_vs_mock'" scripts/run_plot_single_dataset_origvsmock.R Rout/run_plot_single_dataset_origvsmock_$(1)_$(2).Rout
endef
$(foreach k,$(FILT), $(foreach i,$(Dsb),$(eval $(call origvsmockrule_filt,$(i),$(k)))))
$(foreach k,$(FILT), $(foreach i,$(DSbulkb),$(eval $(call origvsmockrule_filt,$(i),$(k)))))
$(foreach k,$(FILT), $(foreach i,$(Dsbsim),$(eval $(call origvsmockrule_filt,$(i),$(k)))))

## ------------------ Summary plots, mostly across mock data sets --------------------- ##
## ------------------------------------------------------------------------------------ ##
$(multidsfigdir)/truefpr/summary_truefpr.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dss),$Y_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='' summarytype='truefpr' dtpext='' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr.Rout

$(multidsfigdir)/pca/summary_pca.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dss),$Y_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='' summarytype='pca' dtpext='' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca.Rout

$(multidsfigdir)/timing/summary_timing.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dss),$Y_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='' summarytype='timing' dtpext='' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing.Rout

$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(Dsb),$Y_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${Dsbc}' filt='' summarytype='crossmethod_consistency' dtpext='' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency.Rout

$(multidsfigdir)/relfprtpr/summary_relfprtpr.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(Dsb),$Y_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${Dsbc}' filt='' summarytype='relfprtpr' dtpext='' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr.Rout

$(multidsfigdir)/fracNA/summary_fracNA.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dss),$Y))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='' summarytype='fracNA' dtpext='' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA.Rout

define summaryrule_truefpr
$(multidsfigdir)/truefpr/summary_truefpr_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dss),$Y_$(1)_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='truefpr' dtpext='' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr,$(k))))

define summaryrule_pca
$(multidsfigdir)/pca/summary_pca_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dss),$Y_$(1)_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='pca' dtpext='' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_pca,$(k))))

define summaryrule_timing
$(multidsfigdir)/timing/summary_timing_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dss),$Y_$(1)_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='timing' dtpext='' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_timing,$(k))))

define summaryrule_crossmethod_consistency
$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency_$(1).rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(Dsb),$Y_$(1)_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${Dsbc}' filt='$(1)' summarytype='crossmethod_consistency' dtpext='' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency,$(k))))

define summaryrule_relfprtpr
$(multidsfigdir)/relfprtpr/summary_relfprtpr_$(1).rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(Dsb),$Y_$(1)_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${Dsbc}' filt='$(1)' summarytype='relfprtpr' dtpext='' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr,$(k))))

define summaryrule_fracna
$(multidsfigdir)/fracNA/summary_fracNA_$(1).rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dss),$Y_$(1)))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dssc}' filt='$(1)' summarytype='fracNA' dtpext='' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_fracna,$(k))))

define summaryrule_filtering
$(multidsfigdir)/filtering/summary_filtering_$(1).rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsb),$Y_$(1)))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsb),$Y))) scripts/run_plot_multi_dataset_summarization.R scripts/summarize_filtering.R \
include_datasets.mk $(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(Dsb),$Y)))
	$R "--args datasets='${Dsbc}' filt='$(1)' summarytype='filtering' dtpext='' figdir='$(multidsfigdir)/filtering' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_filtering_$(1).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering,$(k))))

## ------------------ Summary plots, across mock data sets (bulk) --------------------- ##
## ------------------------------------------------------------------------------------ ##
$(multidsfigdir)/truefpr/summary_truefpr_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dssbulk),$Y_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='' summarytype='truefpr' dtpext='_bulk' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr_bulk.Rout

$(multidsfigdir)/pca/summary_pca_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dssbulk),$Y_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='' summarytype='pca' dtpext='_bulk' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca_bulk.Rout

$(multidsfigdir)/timing/summary_timing_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dssbulk),$Y_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='' summarytype='timing' dtpext='_bulk' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing_bulk.Rout

$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency_bulk.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(DSbulkb),$Y_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${DSbulkb}' filt='' summarytype='crossmethod_consistency' dtpext='_bulk' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency_bulk.Rout

$(multidsfigdir)/relfprtpr/summary_relfprtpr_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(DSbulkb),$Y_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${DSbulkb}' filt='' summarytype='relfprtpr' dtpext='_bulk' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr_bulk.Rout

$(multidsfigdir)/fracNA/summary_fracNA_bulk.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dssbulk),$Y))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='' summarytype='fracNA' dtpext='_bulk' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA_bulk.Rout

define summaryrule_truefpr_bulk
$(multidsfigdir)/truefpr/summary_truefpr_$(1)_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dssbulk),$Y_$(1)_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='$(1)' summarytype='truefpr' dtpext='_bulk' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr_bulk,$(k))))

define summaryrule_pca_bulk
$(multidsfigdir)/pca/summary_pca_$(1)_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dssbulk),$Y_$(1)_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='$(1)' summarytype='pca' dtpext='_bulk' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_pca_bulk,$(k))))

define summaryrule_timing_bulk
$(multidsfigdir)/timing/summary_timing_$(1)_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dssbulk),$Y_$(1)_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='$(1)' summarytype='timing' dtpext='_bulk' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_timing_bulk,$(k))))

define summaryrule_crossmethod_consistency_bulk
$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency_$(1)_bulk.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(DSbulkb),$Y_$(1)_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${DSbulkb}' filt='$(1)' summarytype='crossmethod_consistency' dtpext='_bulk' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency_bulk,$(k))))

define summaryrule_relfprtpr_bulk
$(multidsfigdir)/relfprtpr/summary_relfprtpr_$(1)_bulk.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(DSbulkb),$Y_$(1)_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${DSbulkb}' filt='$(1)' summarytype='relfprtpr' dtpext='_bulk' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr_bulk,$(k))))

define summaryrule_fracna_bulk
$(multidsfigdir)/fracNA/summary_fracNA_$(1)_bulk.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dssbulk),$Y_$(1)))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dssbulk}' filt='$(1)' summarytype='fracNA' dtpext='_bulk' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_fracna_bulk,$(k))))

define summaryrule_filtering_bulk
$(multidsfigdir)/filtering/summary_filtering_$(1)_bulk.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(DSbulkb),$Y_$(1)))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(DSbulkb),$Y))) scripts/run_plot_multi_dataset_summarization.R scripts/summarize_filtering.R \
include_datasets.mk $(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(Dsbulkb),$Y)))
	$R "--args datasets='${DSbulkb}' filt='$(1)' summarytype='filtering' dtpext='_bulk' figdir='$(multidsfigdir)/filtering' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_filtering_$(1)_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering_bulk,$(k))))

## ------------------- Summary plots, across mock data sets (sim) --------------------- ##
## ------------------------------------------------------------------------------------ ##
$(multidsfigdir)/truefpr/summary_truefpr_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dsssim),$Y_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='' summarytype='truefpr' dtpext='_sim' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr_sim.Rout

$(multidsfigdir)/trueperformance/summary_trueperformance_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/performance_realtruth/, $(foreach Y,$(Dsbsim),$Y_performance_realtruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_trueperformance.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='' summarytype='trueperformance' dtpext='_sim' figdir='$(multidsfigdir)/trueperformance' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_trueperformance_sim.Rout

$(multidsfigdir)/pca/summary_pca_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dsssim),$Y_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='' summarytype='pca' dtpext='_sim' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca_sim.Rout

$(multidsfigdir)/timing/summary_timing_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dsssim),$Y_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='' summarytype='timing' dtpext='_sim' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing_sim.Rout

$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency_sim.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(Dsbsim),$Y_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='' summarytype='crossmethod_consistency' dtpext='_sim' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency_sim.Rout

$(multidsfigdir)/relfprtpr/summary_relfprtpr_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(Dsbsim),$Y_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='' summarytype='relfprtpr' dtpext='_sim' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr_sim.Rout

$(multidsfigdir)/fracNA/summary_fracNA_sim.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsssim),$Y))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='' summarytype='fracNA' dtpext='_sim' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA_sim.Rout

define summaryrule_truefpr_sim
$(multidsfigdir)/truefpr/summary_truefpr_$(1)_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(Dsssim),$Y_$(1)_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='$(1)' summarytype='truefpr' dtpext='_sim' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr_sim,$(k))))

define summaryrule_trueperformance_sim
$(multidsfigdir)/trueperformance/summary_trueperformance_$(1)_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/performance_realtruth/, $(foreach Y,$(Dsbsim),$Y_$(1)_performance_realtruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_trueperformance.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='$(1)' summarytype='trueperformance' dtpext='_sim' figdir='$(multidsfigdir)/trueperformance' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_trueperformance_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_trueperformance_sim,$(k))))

define summaryrule_pca_sim
$(multidsfigdir)/pca/summary_pca_$(1)_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(Dsssim),$Y_$(1)_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='$(1)' summarytype='pca' dtpext='_sim' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_pca_sim,$(k))))

define summaryrule_timing_sim
$(multidsfigdir)/timing/summary_timing_$(1)_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(Dsssim),$Y_$(1)_timing))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='$(1)' summarytype='timing' dtpext='_sim' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_timing_sim,$(k))))

define summaryrule_crossmethod_consistency_sim
$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency_$(1)_sim.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(Dsbsim),$Y_$(1)_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='$(1)' summarytype='crossmethod_consistency' dtpext='_sim' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency_sim,$(k))))

define summaryrule_relfprtpr_sim
$(multidsfigdir)/relfprtpr/summary_relfprtpr_$(1)_sim.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(Dsbsim),$Y_$(1)_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='$(1)' summarytype='relfprtpr' dtpext='_sim' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr_sim,$(k))))

define summaryrule_fracna_sim
$(multidsfigdir)/fracNA/summary_fracNA_$(1)_sim.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsssim),$Y_$(1)))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk
	$R "--args datasets='${Dsssimc}' filt='$(1)' summarytype='fracNA' dtpext='_sim' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_fracna_sim,$(k))))

define summaryrule_filtering_sim
$(multidsfigdir)/filtering/summary_filtering_$(1)_sim.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsbsim),$Y_$(1)))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(Dsbsim),$Y))) scripts/run_plot_multi_dataset_summarization.R scripts/summarize_filtering.R \
include_datasets.mk $(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(Dsbsim),$Y)))
	$R "--args datasets='${Dsbsimc}' filt='$(1)' summarytype='filtering' dtpext='_sim' figdir='$(multidsfigdir)/filtering' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_filtering_$(1)_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering_sim,$(k))))

## --------------------------- Summary plots, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk
	$R "--args datasets='${Dsbc}' filt='' dtpext='' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock.Rout

$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_sim.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsbsim),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk
	$R "--args datasets='${Dsbsimc}' filt='' dtpext='_sim' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock_sim.Rout

$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_bulk.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(DSbulkb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk
	$R "--args datasets='${DSbulkb}' filt='' dtpext='_bulk' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock_bulk.Rout

define plotrule_summary_origvsmock
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1).rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk 
	$R "--args datasets='${Dsbc}' filt='$(1)' dtpext='' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmock,$(k))))

define plotrule_summary_origvsmock_sim
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1)_sim.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsbsim),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk 
	$R "--args datasets='${Dsbsimc}' filt='$(1)' dtpext='_sim' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmock_sim,$(k))))

define plotrule_summary_origvsmockbulk
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1)_bulk.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(DSbulkb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk 
	$R "--args datasets='${DSbulkb}' filt='$(1)' dtpext='_bulk' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmockbulk,$(k))))



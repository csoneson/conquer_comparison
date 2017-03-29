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
SUMMARYTYPE := truefpr pca fracNA crossmethod_consistency relfprtpr
SUMMARYTYPE2 := filtering
SUMMARYTYPE3 := trueperformance
SUMMARYTYPE4 := timing

.PHONY: all

## Define the default rule
all: plotds plotind plotorigmock \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_plots.rds \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_sim_plots.rds \
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_bulk_plots.rds \
$(addsuffix _plots.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix _plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _bulk_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)))) \
$(addsuffix _bulk_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _bulk_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE2),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _plots.rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _bulk_plots.rds, $(addprefix $(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_, $(foreach Y,$(FILT),$Y))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach K,$(SUMMARYTYPE3),$(K)/summary_$(K)))) \
$(addsuffix _sim_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE3),$(K)/summary_$(K)_$(Y))))) \
$(addsuffix _all_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE4),$(K)/summary_$(K)_$(Y)))))
#$(addsuffix _bulk_plots.rds, $(addprefix $(multidsfigdir)/, $(foreach Y,$(FILT),$(foreach K,$(SUMMARYTYPE4),$(K)/summary_$(K)_$(Y)))))

## Final plots for publication
finalplots: figures/final_figures/final_plots.rds

## Plot original vs mock comparison
plotorigmock: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(Dsb),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(Dsb),$X_$Y)))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(Dsbsim),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(Dsbsim),$X_$Y)))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach X,$(DSbulksignal),$X))) \
$(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(FILT),$(foreach X,$(DSbulksignal),$X_$Y))))

## Plot individual data set results
plotind: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach X,$(DS),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE),$(foreach Y,$(FILT),$(foreach X,$(DS),$k/$X_$Y_$k))))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE4),$(foreach X,$(DSsimsignal),$k/$X_$k)))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/, $(foreach k,$(PLOTTYPE4),$(foreach Y,$(FILT),$(foreach X,$(DSsimsignal),$k/$X_$Y_$k))))) \
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
$(addsuffix _performance.rds, $(addprefix $(realperfdir)/, $(foreach k,$(DSsimsignal),$k))) \
$(addsuffix _performance.rds, $(addprefix $(realperfdir)/, $(foreach j,$(FILT),$(foreach k,$(DSsimsignal),$k_$j))))

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
define dgerule
results/$(1)_$(2)$(4).rds: scripts/apply_$(2).R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds data/$(1).rds
	$R "--args config_file='config/$(1).json' demethod='$(2)' filt='$(3)'" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$(2)$(4).Rout
endef
$(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerule,$(i),$(j),,))))
$(foreach j,$(MTbulk), $(foreach i,$(DSbulk),$(eval $(call dgerule,$(i),$(j),,))))
$(foreach k, $(FILT), $(foreach j,$(MT), $(foreach i,$(DS),$(eval $(call dgerule,$(i),$(j),$(k),_$(k))))))
$(foreach k, $(FILT), $(foreach j,$(MTbulk), $(foreach i,$(DSbulk),$(eval $(call dgerule,$(i),$(j),$(k),_$(k))))))

## ------------------ Prepare COBRAData object for evaluation ------------------------- ##
## ------------------------------------------------------------------------------------ ##
define cobrarule
$(cobradir)/$(1)$(3)_cobra.rds: scripts/prepare_cobra_for_evaluation.R scripts/calculate_gene_characteristics.R \
$(addsuffix $(3).rds, $(addprefix results/$(1)_, $(foreach Y,$(4),$Y))) scripts/prepare_mae.R include_methods.mk
	$R "--args demethods='$(5)' dataset='$(1)' config_file='config/$(1).json' filt='$(2)' resdir='results' outdir='$(cobradir)'" scripts/prepare_cobra_for_evaluation.R Rout/prepare_cobra_for_evaluation_$(1)$(3)$(6).Rout
endef
$(foreach X,$(DS),$(eval $(call cobrarule,$(X),,,$(MT),${MTc},)))
$(foreach X,$(DSbulk),$(eval $(call cobrarule,$(X),,,$(MTbulk),${MTCbulk},_bulk)))
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call cobrarule,$(X),$(k),_$(k),$(MT),${MTc},))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call cobrarule,$(X),$(k),_$(k),$(MTbulk),${MTcbulk},_bulk))))

## ----------------------------- Calculate concordances ------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define concrule
$(concordancedir)/$(1)$(3)_concordances.rds: scripts/calculate_concordances.R $(cobradir)/$(1)$(3)_cobra.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(concordancedir)'" scripts/calculate_concordances.R Rout/calculate_concordances_$(1)$(3).Rout
endef
$(foreach X,$(DS),$(eval $(call concrule,$(X),,)))
$(foreach X,$(DSbulk),$(eval $(call concrule,$(X),,)))
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call concrule,$(X),$(k),_$(k)))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call concrule,$(X),$(k),_$(k)))))

## ------------------------ Calculate relative performances --------------------------- ##
## ------------------------------------------------------------------------------------ ##
define relperfrule
$(relperfdir)/$(1)$(3)_relative_performance.rds: scripts/calculate_relative_performance_all_truths.R $(cobradir)/$(1)$(3)_cobra.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(relperfdir)'" scripts/calculate_relative_performance_all_truths.R Rout/calculate_relative_performance_all_truths_$(1)$(3).Rout
endef
$(foreach X,$(DS),$(eval $(call relperfrule,$(X),,)))
$(foreach X,$(DSbulk),$(eval $(call relperfrule,$(X),,)))
$(foreach k,$(FILT),$(foreach X,$(DS),$(eval $(call relperfrule,$(X),$(k),_$(k)))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(eval $(call relperfrule,$(X),$(k),_$(k)))))

## -------------------------- Calculate true performances ----------------------------- ##
## ------------------------------------------------------------------------------------ ##
define trueperfrule
$(realperfdir)/$(1)$(3)_performance.rds: scripts/calculate_performance_realtruth.R $(cobradir)/$(1)$(3)_cobra.rds data/$(1)_truth.rds
	$R "--args dataset='$(1)' filt='$(2)' cobradir='$(cobradir)' outdir='$(realperfdir)'" scripts/calculate_performance_realtruth.R Rout/calculate_performance_realtruth_$(1)$(3).Rout
endef
$(foreach X,$(DSsimsignal),$(eval $(call trueperfrule,$(X),,)))
$(foreach k,$(FILT),$(foreach X,$(DSsimsignal),$(eval $(call trueperfrule_filt,$(X),$(k),_$(k)))))

## --------------------------- Plots for evaluation ----------------------------------- ##
## ------------------------------------------------------------------------------------ ##
define plotrule
$(singledsfigdir)/$(2)/$(1)$(4)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)$(4)_cobra.rds
	$R "--args dataset='$(1)' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)$(4)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule,$(X),$(Y),,))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule,$(X),$(Y),,))))
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule_filt,$(X),$(Y),$(k),_$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE1),$(eval $(call plotrule_filt,$(X),$(Y),$(k),_$(k))))))

define plotrule2
$(singledsfigdir)/$(2)/$(1)$(4)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R $(cobradir)/$(1)$(4)_cobra.rds \
$(concordancedir)/$(1)$(4)_concordances.rds scripts/help_function_crossmethod_concordance.R
	$R "--args dataset='$(1)' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)$(4)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y),,))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y),,))))
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y),$(k),_$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE2),$(eval $(call plotrule2,$(X),$(Y),$(k),_$(k))))))

define plotrule3
$(singledsfigdir)/$(2)/$(1)$(4)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(relperfdir)/$(1)$(4)_relative_performance.rds
	$R "--args dataset='$(1)' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)$(4)_$(2).Rout
endef
$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y),,))))
$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y),,))))
$(foreach k,$(FILT),$(foreach X,$(DS),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y),$(k),_$(k))))))
$(foreach k,$(FILT),$(foreach X,$(DSbulk),$(foreach Y,$(PLOTTYPE3),$(eval $(call plotrule3,$(X),$(Y),$(k),_$(k))))))

define plotrule4
$(singledsfigdir)/$(2)/$(1)$(4)_$(2)_summary_data.rds: scripts/run_plot_single_dataset_evaluation.R scripts/plot_single_dataset_$(2).R scripts/plot_setup.R \
$(realperfdir)/$(1)$(4)_performance.rds
	$R "--args dataset='$(1)' filt='$(3)' plottype='$(2)' cobradir='$(cobradir)' concordancedir='$(concordancedir)' relperfdir='$(relperfdir)' realperfdir='$(realperfdir)' figdir='$(singledsfigdir)'" scripts/run_plot_single_dataset_evaluation.R Rout/run_plot_single_dataset_evaluation_$(1)$(4)_$(2).Rout
endef
$(foreach X,$(DSsimsignal),$(foreach Y,$(PLOTTYPE4),$(eval $(call plotrule4,$(X),$(Y),,))))
$(foreach k,$(FILT),$(foreach X,$(DSsimsignal),$(foreach Y,$(PLOTTYPE4),$(eval $(call plotrule4,$(X),$(Y),$(k),_$(k))))))

## -------------------- Plots for characterization of data set ------------------------ ##
## ------------------------------------------------------------------------------------ ##
define plotrule_characterization
$(dschardir)/$(1)$(3)_dataset_characteristics_summary_data.rds: scripts/run_plot_dataset_characterization.R scripts/prepare_mae.R \
subsets/$(1)_subsets.rds data/$(1).rds scripts/calculate_gene_characteristics.R  scripts/calculate_cell_characteristics.R
	$R "--args dataset='$(1)' config_file='config/$(1).json' filt='$(2)' cell_cycle_file='data/cell_cycle_geneids.rds' figdir='$(dschardir)'" scripts/run_plot_dataset_characterization.R Rout/run_plot_dataset_characterization_$(1)$(3).Rout
endef
$(foreach i,$(DS),$(eval $(call plotrule_characterization,$(i),,)))
$(foreach i,$(DSbulk),$(eval $(call plotrule_characterization,$(i),,)))
$(foreach k,$(FILT), $(foreach i,$(DS),$(eval $(call plotrule_characterization,$(i),$(k),_$(k)))))
$(foreach k,$(FILT), $(foreach i,$(DSbulk),$(eval $(call plotrule_characterization,$(i),$(k),_$(k)))))

## -------------------- Plots for evaluation, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
define origvsmockrule
$(figdir)/orig_vs_mock/$(1)$(3)_orig_vs_mock_summary_data.rds: $(concordancedir)/$(1)$(3)_concordances.rds \
$(concordancedir)/$(1)mock$(3)_concordances.rds scripts/plot_setup.R  \
scripts/run_plot_single_dataset_origvsmock.R scripts/plot_single_dataset_origvsmock.R
	$R "--args dataset='$(1)' filt='$(2)' concordancedir='$(concordancedir)' figdir='$(figdir)/orig_vs_mock'" scripts/run_plot_single_dataset_origvsmock.R Rout/run_plot_single_dataset_origvsmock_$(1)$(3).Rout
endef
$(foreach i,$(Dsb),$(eval $(call origvsmockrule,$(i),,)))
$(foreach i,$(DSbulksignal),$(eval $(call origvsmockrule,$(i),,)))
$(foreach i,$(Dsbsim),$(eval $(call origvsmockrule,$(i),,)))
$(foreach k,$(FILT), $(foreach i,$(Dsb),$(eval $(call origvsmockrule,$(i),$(k),_$(k)))))
$(foreach k,$(FILT), $(foreach i,$(DSbulksignal),$(eval $(call origvsmockrule,$(i),$(k),_$(k)))))
$(foreach k,$(FILT), $(foreach i,$(Dsbsim),$(eval $(call origvsmockrule,$(i),$(k),_$(k)))))

## ------------------------ Summary plots, across data sets --------------------------- ##
## ------------------------------------------------------------------------------------ ##
define summaryrule_timing
$(multidsfigdir)/timing/summary_timing_plots.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach Y,$(1),$Y_timing))) \
$(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/timing/, $(foreach K,$(FILT),$(foreach Y,$(1),$Y_$(K)_timing)))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_timing.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(2)' filt='$(3)' summarytype='timing' dtpext='' figdir='$(multidsfigdir)/timing' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_timing.Rout
endef
$(eval $(call summaryrule_timing,$(DS),${DSc},${FILTc}))

define summaryrule_fracNA
$(multidsfigdir)/fracNA/summary_fracNA$(1)_plots.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(2),$Y))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach K,$(FILT),$(foreach Y,$(2),$Y_$(K))))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_fracNA.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(4)' summarytype='fracNA' dtpext='$(1)' figdir='$(multidsfigdir)/fracNA' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_fracNA$(1).Rout
endef
$(eval $(call summaryrule_fracNA,_real,$(DSreal),${DSrealc},${FILTc}))
$(eval $(call summaryrule_fracNA,_bulk,$(DSbulk),${DSbulkc},${FILTc}))



define summaryrule_truefpr
$(multidsfigdir)/truefpr/summary_truefpr$(4)$(1)_plots.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/truefpr/, $(foreach Y,$(2),$Y$(4)_truefpr))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_truefpr.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(5)' summarytype='truefpr' dtpext='$(1)' figdir='$(multidsfigdir)/truefpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_truefpr$(4)$(1).Rout
endef
$(eval $(call summaryrule_truefpr,,$(DSrealmock),${DSrealmockc},,))
$(eval $(call summaryrule_truefpr,_sim,$(DSsimmock),${DSsimmockc},,))
$(eval $(call summaryrule_truefpr,_bulk,$(DSbulkmock),${DSbulkmock},,))
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr,,$(DSrealmock),${DSrealmockc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr,_sim,$(DSsimmock),${DSsimmockc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_truefpr,_bulk,$(DSbulkmock),${DSbulkmock},_$(k),$(k))))

define summaryrule_pca
$(multidsfigdir)/pca/summary_pca$(4)$(1)_plots.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_characterization/, $(foreach Y,$(2),$Y$(4)_results_characterization))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_pca.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(5)' summarytype='pca' dtpext='$(1)' figdir='$(multidsfigdir)/pca' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_pca$(4)$(1).Rout
endef
$(eval $(call summaryrule_pca,,$(DSrealmock),${DSrealmockc},,))
$(eval $(call summaryrule_pca,_sim,$(DSsimmock),${DSsimmockc},,))
$(eval $(call summaryrule_pca,_bulk,$(DSbulkmock),${DSbulkmock},,))
$(foreach k,$(FILT),$(eval $(call summaryrule_pca,,$(DSrealmock),${DSrealmockc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_pca,_sim,$(DSsimmock),${DSsimmockc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_pca,_bulk,$(DSbulkmock),${DSbulkmock},_$(k),$(k))))

define summaryrule_crossmethod_consistency
$(multidsfigdir)/crossmethod_consistency/summary_crossmethod_consistency$(4)$(1)_plots.rds: $(addsuffix .rds, $(addprefix $(concordancedir)/, $(foreach Y,$(2),$Y$(4)_concordances))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_crossmethod_consistency.R scripts/help_function_crossmethod_concordance.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(5)' summarytype='crossmethod_consistency' dtpext='$(1)' figdir='$(multidsfigdir)/crossmethod_consistency' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_crossmethod_consistency$(4)$(1).Rout
endef
$(eval $(call summaryrule_crossmethod_consistency,,$(Dsb),${Dsbc},,))
$(eval $(call summaryrule_crossmethod_consistency,_sim,$(Dsbsim),${Dsbsimc},,))
$(eval $(call summaryrule_crossmethod_consistency,_bulk,$(DSbulksignal),${DSbulksignal},,))
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency,,$(Dsb),${Dsbc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency,_sim,$(Dsbsim),${Dsbsimc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_crossmethod_consistency,_bulk,$(DSbulksignal),${DSbulksignal},_$(k),$(k))))

define summaryrule_relfprtpr
$(multidsfigdir)/relfprtpr/summary_relfprtpr$(4)$(1)_plots.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/results_relativetruth/, $(foreach Y,$(2),$Y$(4)_results_relativetruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_relfprtpr.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(5)' summarytype='relfprtpr' dtpext='$(1)' figdir='$(multidsfigdir)/relfprtpr' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_relfprtpr$(4)$(1).Rout
endef
$(eval $(call summaryrule_relfprtpr,,$(Dsb),${Dsbc},,))
$(eval $(call summaryrule_relfprtpr,_sim,$(Dsbsim),${Dsbsimc},,))
$(eval $(call summaryrule_relfprtpr,_bulk,$(DSbulksignal),${DSbulksignal},,))
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr,,$(Dsb),${Dsbc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr,_sim,$(Dsbsim),${Dsbsimc},_$(k),$(k))))
$(foreach k,$(FILT),$(eval $(call summaryrule_relfprtpr,_bulk,$(DSbulksignal),${DSbulksignal},_$(k),$(k))))

define summaryrule_filtering
$(multidsfigdir)/filtering/summary_filtering_$(1)$(2)_plots.rds: $(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(3),$Y_$(1)))) \
$(addsuffix _cobra.rds, $(addprefix $(cobradir)/, $(foreach Y,$(3),$Y))) scripts/run_plot_multi_dataset_summarization.R scripts/summarize_filtering.R \
include_datasets.mk $(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(3),$Y))) scripts/plot_setup.R
	$R "--args datasets='$(4)' filt='$(1)' summarytype='filtering' dtpext='$(2)' figdir='$(multidsfigdir)/filtering' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_filtering_$(1)$(2).Rout
endef
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering,$(k),,$(Dsb),${Dsbc})))
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering,$(k),_bulk,$(DSbulksignal),${DSbulksignal})))
$(foreach k,$(FILT),$(eval $(call summaryrule_filtering,$(k),_sim,$(Dsbsim),${Dsbsimc})))

define summaryrule_trueperformance
$(multidsfigdir)/trueperformance/summary_trueperformance$(4)$(1)_plots.rds: $(addsuffix _summary_data.rds, $(addprefix $(singledsfigdir)/performance_realtruth/, $(foreach Y,$(2),$Y$(4)_performance_realtruth))) \
scripts/run_plot_multi_dataset_summarization.R scripts/summarize_trueperformance.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='$(3)' filt='$(5)' summarytype='trueperformance' dtpext='$(1)' figdir='$(multidsfigdir)/trueperformance' singledsfigdir='$(singledsfigdir)' cobradir='$(cobradir)' dschardir='$(dschardir)' concordancedir='$(concordancedir)'" scripts/run_plot_multi_dataset_summarization.R Rout/run_plot_multi_dataset_summarization_trueperformance$(4)$(1).Rout
endef
$(eval $(call summaryrule_trueperformance,_sim,$(Dsbsim),${Dsbsimc},,))
$(foreach k,$(FILT),$(eval $(call summaryrule_trueperformance,_sim,$(Dsbsim),${Dsbsimc},_$(k),$(k))))

## --------------------------- Summary plots, orig vs mock ---------------------------- ##
## ------------------------------------------------------------------------------------ ##
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_plots.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='${Dsbc}' filt='' dtpext='' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock.Rout

$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_sim_plots.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsbsim),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='${Dsbsimc}' filt='' dtpext='_sim' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock_sim.Rout

$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_bulk_plots.rds: $(addsuffix _orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(DSbulksignal),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='${DSbulksignal}' filt='' dtpext='_bulk' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/run_plot_multi_dataset_origvsmock_bulk.Rout

define plotrule_summary_origvsmock
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1)_plots.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsb),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R
	$R "--args datasets='${Dsbc}' filt='$(1)' dtpext='' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmock,$(k))))

define plotrule_summary_origvsmock_sim
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1)_sim_plots.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(Dsbsim),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R 
	$R "--args datasets='${Dsbsimc}' filt='$(1)' dtpext='_sim' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock_sim.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmock_sim,$(k))))

define plotrule_summary_origvsmockbulk
$(multidsfigdir)/orig_vs_mock/summary_orig_vs_mock_$(1)_bulk_plots.rds: $(addsuffix _$(1)_orig_vs_mock_summary_data.rds, $(addprefix $(figdir)/orig_vs_mock/, $(foreach Y,$(DSbulksignal),$Y))) \
scripts/run_plot_multi_dataset_origvsmock.R include_datasets.mk scripts/plot_setup.R 
	$R "--args datasets='${DSbulksignal}' filt='$(1)' dtpext='_bulk' origvsmockdir='$(figdir)/orig_vs_mock' figdir='$(multidsfigdir)/orig_vs_mock'" scripts/run_plot_multi_dataset_origvsmock.R Rout/plot_summarize_$(1)_orig_vs_mock_bulk.Rout
endef
$(foreach k,$(FILT),$(eval $(call plotrule_summary_origvsmockbulk,$(k))))

## ----------------------------- Final plots ------------------------------------------ ## 
## ------------------------------------------------------------------------------------ ##
figures/final_figures/final_plots.rds: scripts/make_final_plots.R \
$(multidsfigdir)/fracNA/summary_fracNA_plots.rds $(multidsfigdir)/fracNA/summary_fracNA_TPM_1_25p_plots.rds \
$(multidsfigdir)/fracNA/summary_fracNA_bulk_plots.rds $(multidsfigdir)/fracNA/summary_fracNA_TPM_1_25p_bulk_plots.rds \
$(multidsfigdir)/truefpr/summary_truefpr_plots.rds $(multidsfigdir)/truefpr/summary_truefpr_TPM_1_25p_plots.rds \
$(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(DS),$Y))) \
$(addsuffix _dataset_characteristics_summary_data.rds, $(addprefix $(dschardir)/, $(foreach Y,$(DSbulk),$Y)))
	$R "--args figdir='$(figdir)' outdir='$(figdir)/final_figures'" scripts/make_final_plots.R Rout/make_final_plots.Rout

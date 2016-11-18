## Define the version of R and the path to the library
R = R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define the data sets to run the analysis for
DS := GSE45719 GSE41265

## Define the methods to apply
MT := DESeq2 edgeRLRT Wilcoxon

all: $(addsuffix .rds, $(addprefix results/GSE45719_, ${MT})) $(addsuffix .rds, $(addprefix results/GSE41265_, ${MT}))

.SECONDARY:

config/%.json: scripts/generate_config_%.R
	$R scripts/generate_config_$*.R Rout/generate_config_$*.Rout

subsets/%_subsets.rds: data/%.rds config/%.json scripts/generate_subsets.R
	$R "--args config_file='config/$*.json'" scripts/generate_subsets.R Rout/generate_subsets_$*.Rout

#results/GSE45719_%.rds: scripts/apply_%.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
#	$R "--args config_file='config/GSE45719.json' demethod='$*'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_$*.Rout

#results/GSE41265_%.rds: scripts/apply_%.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE41265_subsets.rds
#	$R "--args config_file='config/GSE41265.json' demethod='$*'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE41265_$*.Rout


define myrule
results/$(1)_%.rds: scripts/apply_%.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/$(1)_subsets.rds
	$R "--args config_file='config/$(1).json' demethod='$$*'" scripts/run_diffexpression.R Rout/run_diffexpression_$(1)_$$*.Rout
endef

$(foreach i,$(DS),$(eval $(call myrule,$(i))))


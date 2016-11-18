## Define the version of R and the path to the library
R = R_LIBS=/home/Shared/Rlib/release-3.4-lib/ /usr/local/R/R-3.3.1/bin/R CMD BATCH --no-restore --no-save

## Define the data sets to run the analysis for
DS = GSE45719

## Define the methods to apply
MT = DESeq2 edgeRLRT Wilcoxon

all: results/GSE45719_edgeRLRT.rds results/GSE45719_DESeq2.rds

config/GSE45719.json: scripts/generate_config_GSE45719.R
	$R scripts/generate_config_GSE45719.R Rout/generate_config_GSE45719.Rout

subsets/%_subsets.rds: data/%.rds config/%.json scripts/generate_subsets.R
	$R "--args config_file='config/$*.json'" scripts/generate_subsets.R Rout/generate_subsets_$*.Rout

results/GSE45719_%.rds: scripts/apply_%.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
	$R "--args config_file='config/GSE45719.json' demethod='$*'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_$*.Rout

#$(word 2,$^)
#subsets/GSE45719_subsets.rds: data/GSE45719.rds config/GSE45719.json scripts/generate_subsets.R
#	$R "--args config_file='config/GSE45719.json'" scripts/generate_subsets.R Rout/generate_subsets_GSE45719.Rout

# results/GSE45719_DESeq2.rds: scripts/apply_DESeq2.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='DESeq2'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_DESeq2.Rout
# 
# results/GSE45719_edgeRLRT.rds: scripts/apply_edgeRLRT.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='edgeRLRT'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_edgeRLRT.Rout
# 
# results/GSE45719_Wilcoxon.rds: scripts/apply_Wilcoxon.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='Wilcoxon'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_Wilcoxon.Rout
# 
# results/GSE45719_edgeRQLF.rds: scripts/apply_edgeRQLF.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='edgeRQLF'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_edgeRQLF.Rout
# 
# results/GSE45719_edgeRLRTrobust.rds: scripts/apply_edgeRLRTrobust.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='edgeRLRTrobust'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_edgeRLRTrobust.Rout
# 
# results/GSE45719_edgeRLRTdeconv.rds: scripts/apply_edgeRLRTdeconv.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='edgeRLRTdeconv'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_edgeRLRTdeconv.Rout
# 
# results/GSE45719_MASTcounts.rds: scripts/apply_MASTcounts.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='MASTcounts'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_MASTcounts.Rout
# 
# results/GSE45719_MASTtpm.rds: scripts/apply_MASTtpm.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='MASTtpm'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_MASTtpm.Rout
# 
# results/GSE45719_MASTcountsDetRate.rds: scripts/apply_MASTcountsDetRate.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='MASTcountsDetRate'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_MASTcountsDetRate.Rout
# 
# results/GSE45719_BPSC.rds: scripts/apply_BPSC.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='BPSC'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_BPSC.Rout
# 
# results/GSE45719_monocle.rds: scripts/apply_monocle.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='monocle'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_monocle.Rout
# 
# results/GSE45719_NODES.rds: scripts/apply_NODES.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='NODES'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_NODES.Rout
# 
# results/GSE45719_SAMseq.rds: scripts/apply_SAMseq.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='SAMseq'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_SAMseq.Rout
# 
# results/GSE45719_SCDE.rds: scripts/apply_SCDE.R scripts/prepare_mae.R scripts/run_diffexpression.R subsets/GSE45719_subsets.rds
# 	$R "--args config_file='config/GSE45719.json' demethod='SCDE'" scripts/run_diffexpression.R Rout/run_diffexpression_GSE45719_SCDE.Rout		
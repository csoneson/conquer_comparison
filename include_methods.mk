comma := ,
empty :=
space := $(empty) $(empty)

## All data sets
DS := GSE45719 GSE45719mock GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 UsoskinGSE59739 GSE60749-GPL13112 GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock
## Data sets for which we have both original and mock results (to compare consistency)
DSb := GSE45719 GSE74596 EMTAB2805 GSE60749-GPL13112 GSE48968-GPL13112
## Data sets to include in summary of characteristics (only mock)
Dss := GSE74596mock GSE45719mock EMTAB2805mock GSE60749-GPL13112mock GSE48968-GPL13112mock
Dssc := $(subst $(space),$(comma),$(Dss))

## All methods
MT := edgeRLRT SAMseq Wilcoxon zingeR edgeRQLF NODES BPSC DESeq2 edgeRLRTdeconv MASTcounts MASTcountsDetRate MASTtpm SCDE monocle edgeRLRTrobust voomlimma zingeRauto Seurat DESeq2census edgeRLRTcensus
MTc := $(subst $(space),$(comma),$(MT))

## All filterings
FILT := TPM_1_25p
FILTc := $(subst $(space),$(comma),$(FILT))

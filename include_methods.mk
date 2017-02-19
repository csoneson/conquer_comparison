comma := ,
empty :=
space := $(empty) $(empty)

## All data sets
DS := GSE45719 GSE45719mock GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 UsoskinGSE59739 UsoskinGSE59739mock GSE60749-GPL13112 GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock
DSc := $(subst $(space),$(comma),$(DS))
## Data sets for which we have both original and mock results (to compare consistency)
Dsb := GSE45719 GSE74596 EMTAB2805 GSE60749-GPL13112 GSE48968-GPL13112 UsoskinGSE59739
Dsbc := $(subst $(space),$(comma),$(Dsb))

## Data sets to include in summary of characteristics (only mock)
Dss := GSE74596mock GSE45719mock EMTAB2805mock GSE60749-GPL13112mock GSE48968-GPL13112mock UsoskinGSE59739mock
Dssc := $(subst $(space),$(comma),$(Dss))

## All methods
MT := edgeRLRT SAMseq Wilcoxon zingeR edgeRQLF NODES NODESnofilt BPSC DESeq2 DESeq2nofilt edgeRLRTdeconv MASTcounts MASTcountsDetRate MASTtpm SCDE edgeRLRTrobust voomlimma zingeRauto Seurat Seuratnofilt DESeq2census edgeRLRTcensus zingeRautonofilt monoclecensus monocle D3E
MTc := $(subst $(space),$(comma),$(MT))

## All filterings
FILT := TPM_1_25p
FILTc := $(subst $(space),$(comma),$(FILT))
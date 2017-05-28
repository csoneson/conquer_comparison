comma := ,
empty :=
space := $(empty) $(empty)

## All methods
MT := edgeRLRT SAMseq Wilcoxon edgeRQLF NODES NODESnofilt BPSC DESeq2 DESeq2nofilt edgeRLRTdeconv MASTcpm MASTcpmDetRate MASTtpm MASTtpmDetRate SCDE edgeRLRTrobust voomlimma SeuratBimod SeuratBimodnofilt SeuratBimodIsExpr2 SeuratTobit DESeq2census edgeRLRTcensus monoclecensus monocle D3E limmatrend ROTSvoom ROTScpm ROTStpm metagenomeSeq
MTc := $(subst $(space),$(comma),$(MT))

## Methods to apply to bulk RNA-seq data sets
MTbulk := edgeRLRT SAMseq Wilcoxon edgeRQLF NODES NODESnofilt DESeq2 DESeq2nofilt edgeRLRTdeconv edgeRLRTrobust voomlimma limmatrend
MTcbulk := $(subst $(space),$(comma),$(MTbulk))

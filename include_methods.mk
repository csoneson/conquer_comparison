comma := ,
empty :=
space := $(empty) $(empty)

## Methods run with R 3.3
MT1 := edgeRLRT SAMseq Wilcoxon edgeRQLF NODES NODESnofilt BPSC DESeq2 DESeq2nofilt edgeRLRTdeconv \
MASTcpm MASTcpmDetRate MASTtpm MASTtpmDetRate SCDE edgeRLRTrobust voomlimma SeuratBimod \
SeuratBimodnofilt SeuratBimodIsExpr2 SeuratTobit DESeq2census edgeRLRTcensus monoclecensus \
monocle D3E limmatrend ROTSvoom ROTScpm ROTStpm metagenomeSeq ttest monoclecount
## Methods run with R 3.4
MT2 := scDD DEsingle# zingeRedgeR zingeRedgeRnofilt 

## All methods
MT := $(MT1) $(MT2)
MTc := $(subst $(space),$(comma),$(MT))

## Methods to apply to bulk RNA-seq data sets
MTbulk := edgeRLRT SAMseq Wilcoxon edgeRQLF NODES DESeq2 DESeq2nofilt edgeRLRTdeconv \
edgeRLRTrobust voomlimma limmatrend ttest
MTcbulk := $(subst $(space),$(comma),$(MTbulk))

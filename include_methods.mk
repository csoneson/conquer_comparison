comma := ,
empty :=
space := $(empty) $(empty)

# ## Methods run with R 3.3
MT3.3 := edgeRLRT Wilcoxon \
SAMseq edgeRQLF NODES NODESnofilt BPSC DESeq2 DESeq2nofilt edgeRLRTdeconv \
MASTcpm MASTcpmDetRate MASTtpm MASTtpmDetRate SCDE edgeRLRTrobust voomlimma SeuratBimod \
SeuratBimodnofilt SeuratBimodIsExpr2 SeuratTobit DESeq2census edgeRLRTcensus monoclecensus \
monocle D3E limmatrend ROTSvoom ROTScpm ROTStpm metagenomeSeq ttest monoclecount \
edgeRQLFDetRate DESeq2betapFALSE
 
## All methods
MT := $(MT3.3)# $(MT3.4) $(MT3.4bc3.6)
MTc := $(subst $(space),$(comma),$(MT))

## Methods to apply to bulk RNA-seq data sets
MTbulk := edgeRLRT SAMseq Wilcoxon edgeRQLF NODES DESeq2 DESeq2nofilt edgeRLRTdeconv \
edgeRLRTrobust voomlimma limmatrend ttest
MTcbulk := $(subst $(space),$(comma),$(MTbulk))

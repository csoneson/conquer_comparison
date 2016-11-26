DS := GSE45719

MT := edgeRLRT SAMseq Wilcoxon edgeRZILRT edgeRQLF NODES monocle monoclecounts BPSC DESeq2 edgeRLRTdeconv edgeRLRTrobust MASTcounts MASTcountsDetRate MASTtpm SCDE 
comma := ,
empty :=
space := $(empty) $(empty)
MTc := $(subst $(space),$(comma),$(MT))

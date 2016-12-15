DS := GSE45719 GSE74596 GSE45719mock EMTAB2805 GSE74596mock EMTAB2805mock GSE77847 GSE77847mock GSE63818-GPL16791

MT := edgeRLRT SAMseq Wilcoxon edgeRZILRT edgeRQLF NODES BPSC DESeq2 edgeRLRTdeconv MASTcounts MASTcountsDetRate MASTtpm SCDE monocle edgeRLRTrobust voomlimma
comma := ,
empty :=
space := $(empty) $(empty)
MTc := $(subst $(space),$(comma),$(MT))

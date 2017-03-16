comma := ,
empty :=
space := $(empty) $(empty)

## All data sets
DS := GSE45719 GSE45719mock GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 GSE60749-GPL13112 GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock UsoskinGSE59739 UsoskinGSE59739mock GSE45719sim123 GSE45719sim123mock GSE74596sim123 GSE74596sim123mock GSE48968-GPL13112sim123 GSE48968-GPL13112sim123mock
DSc := $(subst $(space),$(comma),$(DS))
## Data sets for which we have both original and mock results (to compare consistency)
Dsb := GSE45719 GSE74596 EMTAB2805 GSE60749-GPL13112 GSE48968-GPL13112 UsoskinGSE59739 GSE45719sim123 GSE74596sim123 GSE48968-GPL13112sim123
Dsbc := $(subst $(space),$(comma),$(Dsb))

## Data sets to include in summary of characteristics (only mock)
Dss := GSE74596mock GSE45719mock EMTAB2805mock GSE60749-GPL13112mock GSE48968-GPL13112mock UsoskinGSE59739mock
Dssc := $(subst $(space),$(comma),$(Dss))

## Bulk RNA-seq data sets
DSbulk := EGEUV1 EGEUV1mock
DSbulkb := EGEUV1
Dssbulk := EGEUV1mock

## Data sets to simulate from
DSforsim := GSE45719 GSE74596 GSE48968-GPL13112

## Simulated data sets
DSsim := GSE45719sim123 GSE74596sim123 GSE48968-GPL13112sim123
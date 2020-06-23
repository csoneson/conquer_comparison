comma := ,
empty :=
space := $(empty) $(empty)

## All single-cell data sets
DS := GSE45719 GSE45719mock GSE45719sim123 GSE45719sim123mock \
GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 GSE60749-GPL13112 \
GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock UsoskinGSE59739 UsoskinGSE59739mock  \
GSE74596sim123 GSE74596sim123mock GSE60749-GPL13112sim123 GSE60749-GPL13112sim123mock \
GSE62270-GPL17021 GSE62270-GPL17021mock 10XMonoCytoT 10XMonoCytoTmock
DSc := $(subst $(space),$(comma),$(DS))

## Real data sets for which we have both original and mock results (to compare consistency)
Dsb := GSE45719 GSE74596 EMTAB2805 GSE60749-GPL13112 GSE48968-GPL13112 UsoskinGSE59739 GSE62270-GPL17021 10XMonoCytoT
Dsbc := $(subst $(space),$(comma),$(Dsb))

## Simulated data sets for which we have both original and mock results
Dsbsim := GSE45719sim123 GSE74596sim123 GSE60749-GPL13112sim123
Dsbsimc := $(subst $(space),$(comma),$(Dsbsim))

## Data sets to simulate from
DSforsim := GSE45719 GSE74596 GSE48968-GPL13112 GSE60749-GPL13112

## Real signal data sets
DSrealsignal := GSE45719 GSE74596 EMTAB2805 GSE63818-GPL16791 GSE60749-GPL13112 GSE48968-GPL13112 \
UsoskinGSE59739 GSE62270-GPL17021 10XMonoCytoT
DSrealsignalc := $(subst $(space),$(comma),$(DSrealsignal))
## Real mock data sets
DSrealmock := GSE45719mock GSE74596mock EMTAB2805mock GSE60749-GPL13112mock GSE48968-GPL13112mock \
UsoskinGSE59739mock GSE62270-GPL17021mock 10XMonoCytoTmock
DSrealmockc := $(subst $(space),$(comma),$(DSrealmock))
## All real data sets
DSreal := GSE45719 GSE74596 EMTAB2805 GSE63818-GPL16791 GSE60749-GPL13112 GSE48968-GPL13112 UsoskinGSE59739 GSE45719mock \
GSE74596mock EMTAB2805mock GSE60749-GPL13112mock GSE48968-GPL13112mock UsoskinGSE59739mock \
GSE62270-GPL17021 GSE62270-GPL17021mock 10XMonoCytoT 10XMonoCytoTmock
DSrealc := $(subst $(space),$(comma),$(DSreal))

## Simulated signal data sets
DSsimsignal := GSE45719sim123 GSE74596sim123 GSE60749-GPL13112sim123
DSsimsignalc := $(subst $(space),$(comma),$(DSsimsignal))
## Simulated mock data sets
DSsimmock := GSE45719sim123mock GSE74596sim123mock GSE60749-GPL13112sim123mock
DSsimmockc := $(subst $(space),$(comma),$(DSsimmock))
## All simulated data sets
DSsim := GSE45719sim123 GSE45719sim123mock GSE74596sim123 GSE60749-GPL13112sim123 GSE74596sim123mock GSE60749-GPL13112sim123mock
DSsimc := $(subst $(space),$(comma),$(DSsim))

## Bulk signal data sets
DSbulksignal := EGEUV1
DSbulksignalc := $(subst $(space),$(comma),$(DSbulksignal))
## Bulk mock data sets
DSbulkmock := EGEUV1mock
DSbulkmockc := $(subst $(space),$(comma),$(DSbulkmock))
## All bulk data sets
DSbulk := EGEUV1 EGEUV1mock 
DSbulkc := $(subst $(space),$(comma),$(DSbulk))

## Data sets for tSNE plot
DStsne := GSE45719 GSE74596 GSE48968-GPL13112 EMTAB2805 UsoskinGSE59739 GSE63818-GPL16791 \
GSE60749-GPL13112 GSE62270-GPL17021 10XMonoCytoT EGEUV1
DStsnec := $(subst $(space),$(comma),$(DStsne))

comma := ,
empty :=
space := $(empty) $(empty)

## All single-cell data sets
DS := GSE45719 GSE45719mock GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 GSE60749-GPL13112 \
GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock UsoskinGSE59739 UsoskinGSE59739mock GSE45719sim123 \
GSE45719sim123mock GSE74596sim123 GSE74596sim123mock GSE60749-GPL13112sim123 GSE60749-GPL13112sim123mock \
GSE62270-GPL17021 GSE62270-GPL17021mock 10XMonoCytoT 10XMonoCytoTmock# \
#GSE74596scimpute GSE74596scimputemock GSE74596sim123scimpute GSE74596sim123scimputemock \
#GSE74596drimpute GSE74596drimputemock GSE74596sim123drimpute GSE74596sim123drimputemock \
#GSE74596knnsmooth GSE74596knnsmoothmock GSE74596sim123knnsmooth GSE74596sim123knnsmoothmock
DSc := $(subst $(space),$(comma),$(DS))

## All non-imputed single-cell data sets
DSnonimpute := GSE45719 GSE45719mock GSE74596 GSE74596mock EMTAB2805 EMTAB2805mock GSE63818-GPL16791 GSE60749-GPL13112 \
GSE60749-GPL13112mock GSE48968-GPL13112 GSE48968-GPL13112mock UsoskinGSE59739 UsoskinGSE59739mock GSE45719sim123 \
GSE45719sim123mock GSE74596sim123 GSE74596sim123mock GSE60749-GPL13112sim123 GSE60749-GPL13112sim123mock \
GSE62270-GPL17021 GSE62270-GPL17021mock 10XMonoCytoT 10XMonoCytoTmock
DSnonimputec := $(subst $(space),$(comma),$(DSnonimpute))

## scimpute data sets
#DSscimpute := GSE74596scimpute GSE74596scimputemock GSE74596sim123scimpute GSE74596sim123scimputemock
#DSscimputec := $(subst $(space),$(comma),$(DSscimpute))

## drimpute data sets
#DSdrimpute := GSE74596drimpute GSE74596drimputemock GSE74596sim123drimpute GSE74596sim123drimputemock
#DSdrimputec := $(subst $(space),$(comma),$(DSdrimpute))

## knnsmooth data sets
#DSknnsmooth := GSE74596knnsmooth GSE74596knnsmoothmock GSE74596sim123knnsmooth GSE74596sim123knnsmoothmock
#DSknnsmoothc := $(subst $(space),$(comma),$(DSknnsmooth))

## Real data sets for which we have both original and mock results (to compare consistency)
Dsb := GSE45719 GSE74596 EMTAB2805 GSE60749-GPL13112 GSE48968-GPL13112 UsoskinGSE59739 GSE62270-GPL17021 10XMonoCytoT
Dsbc := $(subst $(space),$(comma),$(Dsb))

## Simulated data sets for which we have both original and mock results
Dsbsim := GSE45719sim123 GSE74596sim123 GSE60749-GPL13112sim123
Dsbsimc := $(subst $(space),$(comma),$(Dsbsim))

## Real imputed data sets for which we have both original and mock results
#Dsbscimpute := GSE74596scimpute
#Dsbscimputec := $(subst $(space),$(comma),$(Dsbscimpute))
#Dsbdrimpute := GSE74596drimpute
#Dsbdrimputec := $(subst $(space),$(comma),$(Dsbdrimpute))
#Dsbknnsmooth := GSE74596knnsmooth
#Dsbknnsmoothc := $(subst $(space),$(comma),$(Dsbknnsmooth))

## Simulated imputed data sets for which we have both original and mock results
#Dsbsimscimpute := GSE74596sim123scimpute
#Dsbsimscimputec := $(subst $(space),$(comma),$(Dsbsimscimpute))
#Dsbsimdrimpute := GSE74596sim123drimpute
#Dsbsimdrimputec := $(subst $(space),$(comma),$(Dsbsimdrimpute))
#Dsbsimknnsmooth := GSE74596sim123knnsmooth
#Dsbsimknnsmoothc := $(subst $(space),$(comma),$(Dsbsimknnsmooth))

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

## Real signal imputed data sets
#DSrealsignalscimpute := GSE74596scimpute
#DSrealsignalscimputec := $(subst $(space),$(comma),$(DSrealsignalscimpute))
#DSrealsignaldrimpute := GSE74596drimpute
#DSrealsignaldrimputec := $(subst $(space),$(comma),$(DSrealsignaldrimpute))
#DSrealsignalknnsmooth := GSE74596knnsmooth
#DSrealsignalknnsmoothc := $(subst $(space),$(comma),$(DSrealsignalknnsmooth))
## Real mock imputed data sets
#DSrealmockscimpute := GSE74596scimputemock
#DSrealmockscimputec := $(subst $(space),$(comma),$(DSrealmockscimpute))
#DSrealmockdrimpute := GSE74596drimputemock
#DSrealmockdrimputec := $(subst $(space),$(comma),$(DSrealmockdrimpute))
#DSrealmockknnsmooth := GSE74596knnsmoothmock
#DSrealmockknnsmoothc := $(subst $(space),$(comma),$(DSrealmockknnsmooth))
## All real imputed data sets
#DSrealscimpute := GSE74596scimpute GSE74596scimputemock
#DSrealscimputec := $(subst $(space),$(comma),$(DSrealscimpute))
#DSrealdrimpute := GSE74596drimpute GSE74596drimputemock
#DSrealdrimputec := $(subst $(space),$(comma),$(DSrealdrimpute))
#DSrealknnsmooth := GSE74596knnsmooth GSE74596knnsmoothmock
#DSrealknnsmoothc := $(subst $(space),$(comma),$(DSrealknnsmooth))

## Simulated signal data sets
DSsimsignal := GSE45719sim123 GSE74596sim123 GSE60749-GPL13112sim123
DSsimsignalc := $(subst $(space),$(comma),$(DSsimsignal))
## Simulated mock data sets
DSsimmock := GSE45719sim123mock GSE74596sim123mock GSE60749-GPL13112sim123mock
DSsimmockc := $(subst $(space),$(comma),$(DSsimmock))
## All simulated data sets
DSsim := GSE45719sim123 GSE74596sim123 GSE60749-GPL13112sim123 GSE45719sim123mock GSE74596sim123mock GSE60749-GPL13112sim123mock
DSsimc := $(subst $(space),$(comma),$(DSsim))

## Simulated signal imputed data sets 
#DSsimsignalscimpute := GSE74596sim123scimpute
#DSsimsignalscimputec := $(subst $(space),$(comma),$(DSsimsignalscimpute))
#DSsimsignaldrimpute := GSE74596sim123drimpute
#DSsimsignaldrimputec := $(subst $(space),$(comma),$(DSsimsignaldrimpute))
#DSsimsignalknnsmooth := GSE74596sim123knnsmooth
#DSsimsignalknnsmoothc := $(subst $(space),$(comma),$(DSsimsignalknnsmooth))
## Simulated mock imputed data sets
#DSsimmockscimpute := GSE74596sim123scimputemock
#DSsimmockscimputec := $(subst $(space),$(comma),$(DSsimmockscimpute))
#DSsimmockdrimpute := GSE74596sim123drimputemock
#DSsimmockdrimputec := $(subst $(space),$(comma),$(DSsimmockdrimpute))
#DSsimmockknnsmooth := GSE74596sim123knnsmoothmock
#DSsimmockknnsmoothc := $(subst $(space),$(comma),$(DSsimmockknnsmooth))
## All simulated imputed data sets
#DSsimscimpute := GSE74596sim123scimpute GSE74596sim123scimputemock
#DSsimscimputec := $(subst $(space),$(comma),$(DSsimscimpute))
#DSsimdrimpute := GSE74596sim123drimpute GSE74596sim123drimputemock
#DSsimdrimputec := $(subst $(space),$(comma),$(DSsimdrimpute))
#DSsimknnsmooth := GSE74596sim123knnsmooth GSE74596sim123knnsmoothmock
#DSsimknnsmoothc := $(subst $(space),$(comma),$(DSsimknnsmooth))

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
DStsne := GSE74596 GSE45719 GSE48968-GPL13112 EMTAB2805 UsoskinGSE59739 GSE63818-GPL16791 \
GSE60749-GPL13112 GSE62270-GPL17021 10XMonoCytoT EGEUV1
DStsnec := $(subst $(space),$(comma),$(DStsne))

impute_dropouts <- function(count, tpm, condt, avetxlength, imputationmethod) {
  if (imputationmethod == "scimpute") {
    source("scripts/scimpute_dropouts.R")
    imputed <- scimpute_dropouts(count = count, tpm = tpm, condt = condt, 
                                 avetxlength = avetxlength)
  } else if (imputationmethod == "drimpute") {
    source("scripts/drimpute_dropouts.R")
    imputed <- drimpute_dropouts(count = count, tpm = tpm, condt = condt, 
                                 avetxlength = avetxlength)
  } else if (imputationmethod == "knnsmooth") {
    source("scripts/knnsmooth_dropouts.R")
    imputed <- knnsmooth_dropouts(count = count, tpm = tpm, condt = condt,
                                  avetxlength = avetxlength)
  }
  imputed
}
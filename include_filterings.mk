comma := ,
empty :=
space := $(empty) $(empty)

## All filterings
FILT := TPM_1_25p
FILTc := $(subst $(space),$(comma),$(FILT))


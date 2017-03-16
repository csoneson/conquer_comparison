args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(config_file)

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(survey))

config <- fromJSON(file = config_file)

print(config)

mae <- readRDS(config$mae)
pdata <- Biobase::pData(mae)

groupid <- config$groupid
keepgroups <- config$keepgroups
sizes <- config$sizes
nreps <- config$nreps

if (length(groupid) > 1) {
  pdata[, paste(groupid, collapse = ".")] <- 
    as.character(interaction(as.data.frame(pdata[, groupid])))
  groupid <- paste(groupid, collapse = ".")
}

## Keep only two groups
if (is.null(keepgroups)) 
  keepgroups <- levels(factor(pdata[, groupid]))[1:2]
keepsamples <- rownames(pdata[pdata[, groupid] %in% keepgroups, , drop = FALSE])
pdata <- droplevels(pdata[match(keepsamples, rownames(pdata)), , drop = FALSE])

condt <- as.character(pdata[, groupid])
names(condt) <- rownames(pdata)

ngroups <- nlevels(factor(condt))
message("Considering the following ", ngroups, ifelse(ngroups == 1, " group: ", " groups: "), 
        paste(levels(factor(condt)), collapse = " vs "))

names(sizes) <- sizes
names(nreps) <- sizes
set.seed(config$seed)

## Subset the data set to the largest sample size, to make sure that all 
## smaller data sets are indeed subsets of the largest data set.
if (length(unique(condt)) == 1) {
  condt <- condt[sort(sample(1:length(condt), 2 * max(sizes)))]
} else {
  condt <- condt[sort(stratsample(as.character(condt), 
                                  structure(rep(max(sizes), 2), 
                                            names = levels(factor(condt)))))]
}

keep_tmp <- lapply(sizes, function(sz) {
  unique(t(sapply(1:nreps[as.character(sz)], function(i) {
    if (length(unique(condt)) == 1) {
      tmpn <- names(condt)
      condt2 <- paste0(condt, ".", sample(rep(c("1", "2"), ceiling(length(condt)/2)))[1:length(condt)])
      names(condt2) <- tmpn
      ngroups <- 2
    } else {
      condt2 <- condt
    }
    smp <- names(condt2)[sort(stratsample(as.character(condt2), 
                                          structure(rep(sz, ngroups), 
                                                    names = levels(factor(condt2)))))]
    cdt <- condt2[smp]
    paste(smp, cdt, sep = "___")
  })))
})

keep_samples <- lapply(keep_tmp, function(w) {
  rbind(apply(w, 2, function(s) sapply(strsplit(s, "___"), .subset, 1)))})
out_condition <- lapply(keep_tmp, function(w) {
  rbind(apply(w, 2, function(s) sapply(strsplit(s, "___"), .subset, 2)))})

saveRDS(list(keep_samples = keep_samples, out_condition = out_condition), 
        file = config$subfile)

sessionInfo()

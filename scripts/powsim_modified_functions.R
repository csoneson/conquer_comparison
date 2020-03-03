## Modified functions from the powsim package, in order to correct the edgeR
## size factors in the parameter estimation
estimateNBParamMod <- function(countData, cData = NULL, design = NULL,
                               RNAseq = "singlecell",
                               estFramework = "edgeR",
                               sigma = 1.96) {
  if (RNAseq == "singlecell") {
    if (estFramework == "edgeR") { #use edgeR
      res = .getDist.edgeR.sc.mod(countData, cData, design, sigma)
    }
  }
  res2 <- c(res, list(RNAseq = RNAseq, estFramework = estFramework, sigma = sigma))
  attr(res2, "param.type") <- "estimated"
  return(res2)
}

.getDist.edgeR.sc.mod <- function(countData = countData, cData = cData, 
                                  design = design, sigma = sigma){
  
  # kick out empty samples and lowly expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  DetectG <- rowMeans(countData) > 0
  countData <- countData[DetectG,fullS]
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))
  
  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }
  
  # scran size factors
  sce <- powsim:::.scran.calc(cnts = countData)
  
  # kick out zero size factor samples
  sf <- scater::sizeFactors(sce)
  sce2 <- suppressWarnings(scater::newSCESet(countData = data.frame(countData[, sf > 0])))
  scater::sizeFactors(sce2) <- sf[sf > 0]
  estS <- sum(sf > 0, na.rm = T)
  
  # convert to edgeR object
  dge1 <- powsim:::.convertToedgeR(sce2)
  
  if (is.null(cData)) {
    cData <- NULL
    dge <- dge1
  } else {
    cData <- as.data.frame(cData, stringsAsFactors = FALSE)
    cData <- cData[sf > 0,]
    dge <- edgeR::DGEList(counts = dge1$counts, lib.size = dge1$samples$lib.size,
                          norm.factors = dge1$samples$norm.factors, samples = cData)
  }
  
  libsize <- dge$samples$lib.size
  names(libsize) <- rownames(dge$samples)
  
  # dropout and expressed component
  if (!is.null(design)) {
    modelmat <- stats::model.matrix(design, data = cData)
    dge2 <- edgeR::estimateDisp(y = dge, design = modelmat, verbose = FALSE)
  }
  if (is.null(design)) {
    dge2 <- edgeR::estimateDisp(y = dge, verbose = FALSE)
  }
  ## Modification: change normalization factors to size factors
  sf <- dge2$samples$norm.factors * dge2$samples$lib.size
  sf <- sf/exp(mean(log(sf)))
  mu2 <- rowMeans(dge2$counts/sf)
  ## End modification
  phi.g <- dge2$tagwise.dispersion
  phi.c <- dge2$common.dispersion
  size <- 1/phi.g
  
  estG <- length(mu2)
  # dropout
  nsamples <- ncol(dge2$counts)
  counts0 <- dge2$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  
  # mean, dispersion and size
  lmu2 <- log2(mu2 + 1)
  ldisp <- log2(phi.g)
  lsize <- log2(size)
  
  ## Modification: change normalization factors to size factors
  sf <- dge2$samples$norm.factors * dge2$samples$lib.size
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- rownames(dge$samples)
  ## End modification
  
  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu2, y = p0, constraint = 'decrease', 
                         nknots = 20, print.warn = FALSE, print.mesg = FALSE)
  cobs.sim <- runif(1000, min(lmu2), max(lmu2))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[, "fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[, "fit"])
  cobs.predict <- cobs.predict[cobs.predict[, "fit"] < 0.05, ]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[, "fit"]), "z"])
  
  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu2, nsigma = sigma)
  
  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu2, nsigma = sigma)
  
  list(seqDepth = libsize, means = mu2, dispersion = phi.g, common.dispersion = phi.c, 
       size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, 
       p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, 
       totalG = totalG, estS = estS, estG = estG)
}
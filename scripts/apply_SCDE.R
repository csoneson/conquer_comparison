suppressPackageStartupMessages(library(scde))

run_SCDE <- function(L) {
  message("scde")
  session_info <- sessionInfo()
  timing <- system.time({
    intcount <- apply(L$count, 2, function(x) {storage.mode(x) <- 'integer'; x})
    o.ifm <- scde.error.models(counts = intcount, groups = L$condt, n.cores = 1,
                               threshold.segmentation = TRUE, 
                               save.crossfit.plots = FALSE, save.model.plots = FALSE,
                               verbose = 0)
    valid.cells <- o.ifm$corr.a > 0
    table(valid.cells)
    o.ifm <- o.ifm[valid.cells, ]
    o.prior <- scde.expression.prior(models = o.ifm, counts = intcount[, valid.cells], 
                                     length.out = 400, show.plot = FALSE)
    grp <- factor(L$condt[which(valid.cells)])
    names(grp) <- rownames(o.ifm)
    ediff <- scde.expression.difference(o.ifm, intcount[, valid.cells], o.prior, 
                                        groups = grp, n.randomizations = 100, 
                                        n.cores = 1, verbose = 0)
    p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail = FALSE)
  })
  
  hist(p.values.adj, 50)
  
  list(session_info = session_info,
       timing = timing,
       res = ediff,
       df = data.frame(padj = p.values.adj,
                       score = abs(ediff$Z), 
                       row.names = rownames(ediff)))
}
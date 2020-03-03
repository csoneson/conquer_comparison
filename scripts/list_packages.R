lf <- list.files("Rout")
lf <- lf[-grep("zingeR", lf)]
lf <- lf[-grep("NODESnofilt", lf)]
all_packages <- c()
for (f in lf) {
  x <- readLines(paste0("Rout/", f))
  idx1 <- which(x == "other attached packages:")
  idx2 <- which(x == "loaded via a namespace (and not attached):")
  if (length(idx1) != 0 & length(idx2) != 0) {
    all_packages <- 
      unique(c(all_packages, 
               do.call(c, lapply((idx1 + 1):(idx2 - 2), function(i) {
                 grep("\\[", setdiff(setdiff(strsplit(x[i], " ")[[1]], " "), ""), 
                      value = TRUE, invert = TRUE)
               }))))
  }
}
write.table(sort(all_packages), file = "all_packages.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
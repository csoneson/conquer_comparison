run_D3E <- function(L) {
  message("D3E")
  session_info <- sessionInfo()
  timing <- system.time({
    tmp <- cbind(GeneID = rownames(L$count), L$count)
    colnames(tmp) <- c("GeneID", make.names(L$condt))
    rnb <- paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), "_", round(runif(1) * 1e8))
    message(rnb)
    write.table(tmp, file = paste0("tmp/d3e_", rnb, ".txt"), row.names = FALSE,
                col.names = TRUE, quote = FALSE, sep = "\t")
    cmd <- sprintf("python2 software/D3E/D3ECmd.py %s %s %s %s -m 0 -t 0 -z 0 -n 1 -v",
                   paste0("tmp/d3e_", rnb, ".txt"),
                   paste0("tmp/d3e_", rnb, ".out"),
                   levels(factor(make.names(L$condt)))[1], 
                   levels(factor(make.names(L$condt)))[2])
    message(cmd)
    system(cmd)
    res <- read.delim(paste0("tmp/d3e_", rnb, ".out"), header = TRUE,
                      as.is = TRUE, row.names = 1)
    
  })
  
  hist(res$p.value, 50)
  
  list(session_info = session_info,
       timing = timing,
       res = res,
       df = data.frame(pval = res$p.value,
                       padj = p.adjust(res$p.value, method = "BH"), 
                       row.names = rownames(res)))
}
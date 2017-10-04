plot_runfailure <- function(cobra, runstatuslist, colvec, exts, summary_data = list()) {
  runstatus <- unlist(runstatuslist)
  runstatus <- data.frame(method = names(runstatus), result = runstatus) %>%
    tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(exts, "", method))
  
  ## Extract all method/instance combinations with results from the cobra object
  success <- union(union(colnames(pval(cobra)), colnames(padj(cobra))),
                   colnames(iCOBRA::score(cobra)))
  df <- data.frame(method = success, stringsAsFactors = FALSE) %>%
    tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(exts, "", method)) %>%
    dplyr::mutate(incobra = "yes")
  
  ## Get all instances
  instc <- data.frame(instc = grep("tested", colnames(truth(cobra)), value = TRUE), 
                      stringsAsFactors = FALSE) %>%
    tidyr::separate(instc, into = c("status", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells)))
  instc <- do.call(rbind, lapply(unique(df$method), function(w) {
    instc %>% dplyr::mutate(method = w)
  }))
  
  df <- dplyr::full_join(df, instc, by = c("method", "ncells", "repl"))
  df$incobra[is.na(df$incobra)] <- "no"
  
  runstatus <- dplyr::full_join(runstatus, df, by = c("method", "ncells", "repl"))
  
  ## Summarize failure rates
  runstatus_sum <- runstatus %>% dplyr::group_by(method) %>% 
    dplyr::summarise(failure_rate = sum(result == "failure")/sum(!is.na(result)))
  
  summary_data$runstatus <- runstatus
  summary_data$runstatus_sum <- runstatus_sum
  return(invisible(summary_data))
}


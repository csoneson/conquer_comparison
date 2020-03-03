plot_runfailure <- function(runstatuslist, colvec, exts, summary_data = list()) {
  runstatus <- unlist(runstatuslist)
  runstatus <- data.frame(method = names(runstatus), result = runstatus) %>%
    tidyr::separate(method, into = c("method", "ncells", "repl"), sep = "\\.") %>%
    dplyr::mutate(ncells = as.numeric(as.character(ncells))) %>%
    dplyr::mutate(method = gsub(exts, "", method))

  ## Summarize failure rates
  runstatus_sum <- runstatus %>% dplyr::group_by(method) %>% 
    dplyr::summarise(failure_rate = sum(result == "failure")/sum(!is.na(result)))
  
  summary_data$runstatus <- runstatus
  summary_data$runstatus_sum <- runstatus_sum
  
  return(invisible(summary_data))
}


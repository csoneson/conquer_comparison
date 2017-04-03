suppressPackageStartupMessages(library(testthat))

test_results <- test_dir("unit_tests", reporter = "summary")

test_results

sessionInfo()
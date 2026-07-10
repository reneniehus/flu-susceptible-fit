#!/usr/bin/env Rscript
# Run the playground's test suite from the project root:
#   Rscript run_tests.R
# Tests are fully self-contained: they SIMULATE their own data (no internet, no external files) and
# check the engine's invariants, the truth/observed contract, and that every analysis tool recovers
# the truth it is meant to. See tests/testthat/.
Sys.setenv(PP_ROOT = getwd())
library(testthat)
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

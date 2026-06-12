#!/usr/bin/env Rscript
# Run the project's test suite from the repo root:
#   Rscript run_tests.R
# Tests build the pipeline from the local data/ snapshots (no internet needed) and check the
# data contracts and key invariants. See tests/testthat/.
library(testthat)
library(here)
testthat::test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)

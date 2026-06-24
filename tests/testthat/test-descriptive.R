# Tests for the descriptive curve-feature method (code/01_main_supporting/methods/method_descriptive.R).
# Smooths the slim DK panel offline and checks the common contract, the implied S0, and that the
# descriptive features are well-formed.

sl  <- load_flu_iliplus_slim("DK", path = here::here("data/slim_flu_iliplus.csv"))
fit <- fit_descriptive(sl$ylist, R0 = 1.5, infectious_period_days = 3, smooth_window = 4)
fit$country <- "DK"; fit$seasons <- sl$seasons; fit$season_week <- sl$season_week   # as run_method attaches
K   <- length(sl$ylist)

test_that("descriptive fit returns the common contract (curve only, no mechanistic params)", {
  expect_equal(fit$method, "descriptive")
  expect_equal(fit$convergence, 0L)
  expect_length(fit$mu, K)
  expect_true(all(is.na(fit$params$S0)))                         # no implied S0 -- steepness is the feature
  expect_true(all(is.na(fit$params$c)) && is.na(fit$params$qI))  # no SIR parameters
})

test_that("the smoothed curve tracks the observed data", {
  cors <- vapply(seq_len(K), function(s){
    y <- sl$ylist[[s]]; mu <- fit$mu[[s]]; ok <- is.finite(y) & is.finite(mu)
    if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_
  }, numeric(1))
  expect_gt(median(cors, na.rm = TRUE), 0.9)
})

test_that("the descriptive features are well-formed", {
  smry <- summarise_method_fit(fit)
  expect_equal(nrow(smry), K)
  expect_true(all(smry$auc > 0))
  expect_true(all(smry$peak_height > 0))
  expect_true(all(is.finite(smry$steepness)))
  expect_true(all(is.finite(smry$onset_week)))
})

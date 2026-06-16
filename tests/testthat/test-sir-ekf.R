# Tests for the EKF SIR method (code/01_main_supporting/methods/method_sir_ekf.R). Fits the slim DK
# panel offline and checks the common contract, that the process noise is fitted but stays small
# (so it regularises rather than masks S0), and that the filter reproduces the waves.

sl  <- load_flu_iliplus_slim("DK", path = here::here("data/slim_flu_iliplus.csv"))
fit <- fit_sir_ekf(sl$ylist, R0 = 1.5, infectious_period_days = 3,
                   seed_i0 = 1e-5, n_starts = 2, seed = 1)
K   <- length(sl$ylist)

test_that("EKF fit returns the common contract with plausible parameters", {
  expect_equal(fit$method, "ekf")
  expect_equal(fit$R0, 1.5); expect_equal(fit$seed_i0, 1e-5)     # R0 and seed are FIXED, not fitted
  expect_equal(fit$convergence, 0)
  expect_length(fit$params$S0, K)                                # one susceptibility per season
  expect_length(fit$params$c,  K)                                # one reporting fraction per season
  expect_length(fit$mu, K)
  expect_true(all(fit$params$S0 > 0 & fit$params$S0 < 1))
  expect_gt(median(fit$params$S0), 1 / fit$R0)                   # typical season in the growing regime
  expect_true(all(fit$params$c > 0)); expect_gt(fit$params$phi, 0); expect_gt(fit$params$b, 0)
})

test_that("the process noise is fitted, positive and small (regularises, does not mask S0)", {
  expect_true(is.finite(fit$params$qI) && fit$params$qI > 0)
  expect_lt(fit$params$qI, 1e-2)                                 # stays near its small prior
})

test_that("the EKF reproduces each season's wave", {
  cors <- vapply(seq_len(K), function(s){
    y <- sl$ylist[[s]]; mu <- fit$mu[[s]]; ok <- is.finite(y) & is.finite(mu)
    if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_
  }, numeric(1))
  expect_gt(median(cors, na.rm = TRUE), 0.8)                    # filter tracks, so high correlation
})

test_that("trajectory projects beyond the data window", {
  n  <- length(sl$ylist[[1]])
  expect_length(sir_ekf_trajectory(fit, 1, n_weeks = n + 8), n + 8)
  expect_true(all(sir_ekf_trajectory(fit, 1) >= 0))
})

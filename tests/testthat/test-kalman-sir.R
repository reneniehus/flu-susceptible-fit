# Tests for the EKF-SIR re-implementation (code/01_main_supporting/model_kalman_sir.R).
# Fits flu ILI+ for one country/season offline and checks the fit is sensible.

s   <- kalman_sir_series(models_in, "DK", "2023/2024")
fit <- fit_kalman_sir(s$value, infectious_period_days = 3, n_sub = 7, n_starts = 4, seed = 1)

test_that("EKF-SIR fit converges with epidemiologically plausible parameters", {
  expect_equal(fit$convergence, 0)
  expect_gt(fit$params$R0, 1.0); expect_lt(fit$params$R0, 3.0)   # priors keep R0 in range
  expect_gt(fit$params$S0, 0.2); expect_lt(fit$params$S0, 0.99)
  expect_gt(fit$params$I0, 1e-9)                                  # no degenerate zero seed
})

test_that("EKF-SIR reproduces the observed flu ILI+ wave", {
  ok <- is.finite(s$value) & is.finite(fit$mu_pred)
  expect_gt(cor(s$value[ok], fit$mu_pred[ok]), 0.8)
  # peak timing within a couple of weeks
  peak_obs <- s$season_week[which.max(replace(s$value, !is.finite(s$value), -Inf))]
  peak_fit <- s$season_week[which.max(fit$mu_pred)]
  expect_lt(abs(peak_obs - peak_fit), 6)
})

test_that("deterministic trajectory has the right length and is non-negative", {
  tr <- kalman_sir_trajectory(fit)
  expect_length(tr, length(s$value))
  expect_true(all(tr >= 0))
  expect_length(kalman_sir_trajectory(fit, n_weeks = length(s$value) + 8), length(s$value) + 8) # projects
})

# Tests for the Stan-matched multi-season re-parameterisation of the EKF-SIR
# (code/01_main_supporting/model_kalman_sir.R): R0 fixed, per-season initial conditions,
# shared reporting fraction c and overdispersion phi. Runs offline on the committed slim
# panel (data/slim_flu_iliplus.csv).

sl  <- load_flu_iliplus_slim("DK", path = here::here("data/slim_flu_iliplus.csv"))
fit <- fit_kalman_sir_multiseason(sl$ylist, R0 = 1.3, infectious_period_days = 3,
                                  n_sub = 7, n_starts = 4, seed = 1)
K   <- length(sl$ylist)

test_that("multi-season fit keeps R0 fixed and the parameter shapes match Stan's", {
  expect_equal(fit$convergence, 0)
  expect_equal(fit$R0, 1.3)                                   # R0 is fixed, not fitted
  expect_length(fit$params$S0, K)                            # one initial S per season
  expect_length(fit$params$I0, K)                            # one initial I per season
  expect_length(fit$params$c, 1L)                            # reporting fraction shared
  expect_length(fit$params$phi, 1L)                          # overdispersion shared
})

test_that("per-season initial conditions are epidemiologically plausible", {
  expect_true(all(fit$params$S0 > 0.2 & fit$params$S0 < 0.99))
  expect_true(all(fit$params$I0 > 1e-9))                     # no degenerate zero seed
  expect_gt(fit$params$c, 0)
  expect_gt(fit$params$phi, 0)
})

test_that("the shared-parameter fit reproduces each season's flu ILI+ wave", {
  cors <- vapply(seq_len(K), function(s){
    y <- sl$ylist[[s]]; mu <- fit$fits[[s]]$mu_pred
    ok <- is.finite(y) & is.finite(mu)
    if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_
  }, numeric(1))
  expect_true(median(cors, na.rm = TRUE) > 0.7)             # filtered fit tracks the waves
})

test_that("per-season deterministic trajectory has the right length and is non-negative", {
  tr <- kalman_sir_ms_trajectory(fit, 1)
  expect_length(tr, length(sl$ylist[[1]]))
  expect_true(all(tr >= 0))
  expect_length(kalman_sir_ms_trajectory(fit, 1, n_weeks = length(sl$ylist[[1]]) + 8),
                length(sl$ylist[[1]]) + 8)                  # projects forward
})

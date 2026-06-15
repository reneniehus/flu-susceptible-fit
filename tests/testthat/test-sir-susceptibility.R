# Tests for the deterministic per-season susceptibility fit (fit_sir_susceptibility,
# code/01_main_supporting/model_kalman_sir.R). Fits the committed slim DK panel offline and checks
# the fit is sensible: R0/seed stay fixed, per-season S0 and c are recovered, and the DETERMINISTIC
# SIR (no filter) reproduces the observed waves.

sl  <- load_flu_iliplus_slim("DK", path = here::here("data/slim_flu_iliplus.csv"))
fit <- fit_sir_susceptibility(sl$ylist, R0 = 1.5, infectious_period_days = 3,
                              seed_i0 = 1e-5, n_starts = 3, seed = 1)
K   <- length(sl$seasons)

test_that("susceptibility fit converges with plausible per-season parameters", {
  expect_equal(fit$convergence, 0)
  expect_equal(fit$R0, 1.5); expect_equal(fit$seed_i0, 1e-5)      # R0 and seed are FIXED, not fitted
  expect_length(fit$params$S0, K)                                 # one susceptibility per season
  expect_length(fit$params$c,  K)                                 # one reporting fraction per season
  expect_true(all(fit$params$S0 > 0 & fit$params$S0 < 1))
  expect_gt(median(fit$params$S0), 1 / fit$R0)                    # typical season is in the growing regime
  expect_true(all(fit$params$c > 0))
  expect_gt(fit$params$phi, 0); expect_gt(fit$params$b, 0)        # shared baseline + overdispersion
})

test_that("the deterministic SIR reproduces each season's wave", {
  cors <- vapply(seq_len(K), function(s){
    y <- sl$ylist[[s]]; mu <- fit$mu[[s]]; ok <- is.finite(y) & is.finite(mu)
    if (sum(ok) > 2) cor(y[ok], mu[ok]) else NA_real_
  }, numeric(1))
  expect_gt(median(cors, na.rm = TRUE), 0.75)                    # reconstruction, not filter tracking
})

test_that("trajectory has the right length, is non-negative, and projects forward", {
  n  <- length(sl$ylist[[1]])
  tr <- sir_susceptibility_trajectory(fit, 1)
  expect_length(tr, n)
  expect_true(all(tr >= 0))
  expect_length(sir_susceptibility_trajectory(fit, 1, n_weeks = n + 8), n + 8)  # projects beyond data
})

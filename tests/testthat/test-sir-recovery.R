# Parameter recovery on synthetic truth for the SIR susceptibility methods. The project's central
# claim (documentation/decisions.md) is that with R0, gamma and the seed I0 FIXED, the per-season
# susceptibility S0 is identified by the wave's rise rate r = gamma*(R0*S0 - 1). The real-data
# tests can only check plausibility; here we generate seasons from the generative model itself
# (.sir_season_incidence in sir_core.R) with KNOWN, distinct S0 and check both fitting methods
# get them back -- rank order exactly, absolute values within a calibrated tolerance.
#
# Truth: S0 = 0.75 / 0.80 / 0.85. All must sit in the supercritical regime: with R0 = 1.5 a wave
# only grows when S0 > 1/R0 ~ 0.667, and the peak must land inside the 33-week window for the rise
# rate to be observed (S0 <= 0.70 leaves the peak beyond week 33 and is not identifiable from one
# season). Clean peaks fall at weeks 27/19/15 -- realistic mid-winter timing. Reporting c = 4000
# puts the noisy peaks at ~60-240, the scale of the real panel; baseline b = 2; noise is
# rnbinom(mu, size = phi) whose Var = mu + mu^2/phi is exactly the methods' observation model,
# with phi = 30 (~20% CV at the peak). The noise seed is FIXED, so the test is deterministic.
#
# Measured recovery at this seed (calibration runs, this repo):
#   deterministic: S0_hat - S0_true = +0.0002 / -0.0011 / +0.0003  (max abs 0.0011; tol 0.02)
#   EKF:           S0_hat - S0_true = -0.0071 / -0.0152 / -0.0371  (max abs 0.0371; tol 0.08)
# Across three other noise seeds the deterministic worst case was 0.045 (one inflated-c local mode
# on the slowest season) and the EKF worst case 0.051; rank order was correct in every run. The
# tolerances above cover the fixed seed with a wide margin.

syn_gamma <- 7 / params$susc_infectious_period_days        # per-week recovery rate
syn_beta  <- params$susc_R0 * syn_gamma
S0_true   <- c(0.75, 0.80, 0.85)                           # distinct known truths, increasing
c_true    <- 4000                                          # reporting fraction (same every season)
b_true    <- 2                                             # small off-season baseline
phi_true  <- 30                                            # overdispersion: Var = mu + mu^2/phi
n_wk_syn  <- 33                                            # one season, Aug -> late March

set.seed(20260813)                                         # FIXED noise seed -> deterministic test
ylist_syn <- lapply(S0_true, function(S0){
  inc <- .sir_season_incidence(S0, params$susc_seed_i0, syn_beta, syn_gamma, n_wk_syn, n_sub = 7)
  as.numeric(rnbinom(n_wk_syn, mu = c_true * inc + b_true, size = phi_true))
})
names(ylist_syn) <- sprintf("synthetic_S0_%.2f", S0_true)

det_syn <- fit_sir_deterministic(ylist_syn, R0 = params$susc_R0,
                                 infectious_period_days = params$susc_infectious_period_days,
                                 seed_i0 = params$susc_seed_i0, n_starts = 2, seed = 1)
ekf_syn <- fit_sir_ekf(ylist_syn, R0 = params$susc_R0,
                       infectious_period_days = params$susc_infectious_period_days,
                       seed_i0 = params$susc_seed_i0, n_starts = 2, seed = 1)

test_that("the synthetic truth shows the identification mechanism: higher S0 = steeper, earlier, bigger wave", {
  peaks      <- vapply(ylist_syn, max, numeric(1))
  peak_weeks <- vapply(ylist_syn, which.max, integer(1))
  expect_true(all(diff(peaks) > 0))                        # bigger waves with more susceptibility
  expect_true(all(diff(peak_weeks) < 0))                   # faster rise -> earlier peak
})

test_that("deterministic fit recovers the S0 rank order exactly and each S0 within 0.02", {
  expect_equal(det_syn$convergence, 0)
  expect_equal(order(det_syn$params$S0), order(S0_true))   # relative S0 is the quantity of interest
  expect_lt(max(abs(det_syn$params$S0 - S0_true)), 0.02)   # calibrated: measured max abs error 0.0011
})

test_that("deterministic fit recovers the reporting fraction to within a factor of 2", {
  # c and S0 control different curve features (level vs rise rate), so with the true I0 supplied
  # c must come back near truth -- if S0 errors were being absorbed into c this would fail.
  expect_true(all(det_syn$params$c / c_true > 0.5 & det_syn$params$c / c_true < 2))
})

test_that("EKF recovers the S0 rank order exactly and each S0 within 0.08", {
  # Looser absolute tolerance than the deterministic method: the EKF's process noise and tight
  # regularising priors shrink S0 slightly, but the ranking -- the interpreted quantity -- must hold.
  expect_equal(ekf_syn$convergence, 0)
  expect_equal(order(ekf_syn$params$S0), order(S0_true))
  expect_lt(max(abs(ekf_syn$params$S0 - S0_true)), 0.08)   # calibrated: measured max abs error 0.0371
  expect_lt(ekf_syn$params$qI, 1e-2)                       # process noise regularises, does not absorb the wave
})

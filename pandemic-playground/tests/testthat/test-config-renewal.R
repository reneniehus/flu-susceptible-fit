# Config validation (R/config.R) and the renewal engine (R/renewal.R).

test_that("the default config validates, and validate_config catches broken configs", {
  cfg <- default_config()
  expect_true(validate_config(cfg))
  expect_error(validate_config(within(cfg, {ifr <- 2})))              # IFR out of (0,1)
  bad <- cfg; bad$ascertainment$death_rho <- 1.5
  expect_error(validate_config(bad))
  bad2 <- cfg; bad2$delays$generation_interval <- "not an epidist"
  expect_error(validate_config(bad2))
})

test_that("config_subset keeps every country-indexed field consistent", {
  cfg <- config_subset(default_config(), c("DE", "FR"), n_days = 90)
  expect_equal(cfg$countries, c("DE", "FR"))
  expect_equal(nrow(cfg$geography), 2)
  expect_named(cfg$rt_country, c("DE", "FR"))
  expect_named(cfg$surveillance_quality, c("DE", "FR"))
  expect_equal(cfg$n_days, 90)
  expect_true(validate_config(cfg))
  expect_error(config_subset(cfg, "ZZ"))                             # unknown country
})

test_that("renewal grows when R>1 and dies out when R<1", {
  gi <- discretise(epidist_gamma("gi", 5, 1.7), boundary = "cori")
  seed <- matrix(0, 80, 1); seed[1, 1] <- 20
  up   <- simulate_renewal(80, 1e7, gi, list(step_schedule(0, 2.0)), seed, stochastic = FALSE)
  down <- simulate_renewal(80, 1e7, gi, list(step_schedule(0, 0.6)), seed, stochastic = FALSE)
  expect_gt(up$incidence[60, 1], up$incidence[10, 1])                # growing
  # for R<1 the single seed's generation-interval kernel first ramps (peak ~day 7) then decays; check
  # the sustained decline of successive generations rather than the pre-peak day-2 value
  expect_lt(down$incidence[60, 1], down$incidence[10, 1])            # declining generations
  expect_true(all(diff(down$incidence[10:70, 1]) < 0))              # monotone decay past the seed peak
})

test_that("susceptibles are conserved and depletion turns the epidemic over", {
  gi <- discretise(epidist_gamma("gi", 5, 1.7), boundary = "cori")
  seed <- matrix(0, 200, 1); seed[1, 1] <- 20
  out  <- simulate_renewal(200, 1e5, gi, list(step_schedule(0, 2.5)), seed, stochastic = FALSE)
  expect_true(all(diff(out$susceptible) <= 1e-6))                    # S monotone non-increasing
  expect_true(all(out$susceptible >= 0))
  # cumulative infections + remaining susceptibles ~ N (nothing created or lost)
  expect_equal(sum(out$incidence) + tail(out$susceptible, 1), 1e5, tolerance = 1e-4 * 1e5)
  expect_lt(which.max(out$incidence[, 1]), 200)                      # peaks before the horizon (turnover)
})

test_that("validate_config catches a missing surveillance value and a non-positive dispersion_k", {
  cfg <- default_config()
  bad <- cfg; bad$surveillance_quality <- bad$surveillance_quality[-1]   # drop one country's value
  expect_error(validate_config(bad), "surveillance")                    # the check now actually fires
  bad2 <- cfg; bad2$dispersion_k <- 0
  expect_error(validate_config(bad2), "dispersion")
})

test_that("the Poisson limit (dispersion_k = Inf) is identical to the default draw", {
  cfg <- config_subset(default_config(), c("DE", "FR"), n_days = 100)   # default dispersion_k = Inf
  a <- simulate_pandemic(cfg)
  cfg2 <- cfg; cfg2$dispersion_k <- Inf
  b <- simulate_pandemic(cfg2)
  expect_identical(a$truth$infections, b$truth$infections)
})

test_that("transmission overdispersion (dispersion_k) elevates chain extinction vs Poisson", {
  gi <- discretise(epidist_gamma("gi", 5, 1.7), boundary = "cori")
  establish <- function(k, reps = 80) mean(vapply(seq_len(reps), function(i) {
    set.seed(i); seed <- matrix(0, 100, 1); seed[1, 1] <- 1
    sum(simulate_renewal(100, 1e7, gi, list(step_schedule(0, 2.4)), seed, dispersion = k)$incidence) > 50
  }, logical(1)))
  p_poisson <- establish(Inf); p_super <- establish(0.16)
  expect_gt(p_poisson, 0.7)                            # Poisson: most single-seed chains establish
  expect_lt(p_super, 0.5)                              # superspreading (SARS k): most fizzle out
  expect_gt(p_poisson - p_super, 0.25)                # a large, real difference
})

test_that("exogenous seeding plants infections on the given day", {
  gi <- discretise(epidist_gamma("gi", 5, 1.7), boundary = "cori")
  seed <- matrix(0, 60, 2); seed[1, 1] <- 10; seed[30, 2] <- 5      # strain 2 introduced on day 29
  out  <- simulate_renewal(60, 1e7, gi, list(step_schedule(0, 1.5), step_schedule(0, 2.0)), seed, stochastic = FALSE)
  expect_equal(out$incidence[1, 1], 10)                             # day-0 seed present
  expect_equal(sum(out$incidence[1:30, 2]), 5)                      # strain 2 only from its intro day
  expect_gt(out$incidence[50, 2], 0)                                # and then it grows
})

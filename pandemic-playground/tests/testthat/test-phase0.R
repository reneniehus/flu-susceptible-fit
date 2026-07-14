# Phase 0 tools recover the truth (analysis/phase0_*.R).

test_that("growth->R recovers the true early R on the true infection curve", {
  gi <- discretise(test_input$delays$generation_interval, boundary = "cori")
  inf <- truth_infections(test_sim, "X"); d <- inf[inf$day %in% 20:45, ]
  R <- r_to_R(estimate_growth_rate(d$infections, d$day)$r, gi)
  expect_equal(R, 2.6, tolerance = 0.15)                            # true source R = 2.6
})

test_that("r_to_R inverts the Euler-Lotka relation", {
  gi <- discretise(test_input$delays$generation_interval, boundary = "cori")
  a <- 0:(length(gi) - 1)
  r_star <- uniroot(function(r) 2.0 * sum(gi * exp(-r * a)) - 1, c(1e-4, 1))$root
  expect_equal(r_to_R(r_star, gi), 2.0, tolerance = 1e-3)
})

test_that("growth doubling-time CI is NA when the growth-rate CI straddles zero", {
  set.seed(1); counts <- stats::rpois(30, 20)          # a ~flat series -> r near 0, CI includes 0
  g <- estimate_growth_rate(counts, 0:29)
  if (g$r_lower < 0 && g$r_upper > 0) {
    expect_true(is.na(g$doubling_lower))               # a doubling time is undefined when r may be <= 0
    expect_true(is.na(g$doubling_upper))
  }
  expect_true(g$char_kind %in% c("doubling", "halving"))
})

test_that("the chain-size PMF is a valid distribution (NB and Borel limits)", {
  psum <- function(R, k) sum(exp(vapply(1:2000, .chain_logpmf, numeric(1), R = R, k = k)))
  expect_equal(psum(0.9, 0.5), 1, tolerance = 1e-3)    # negative-binomial offspring
  expect_equal(psum(0.7, Inf), 1, tolerance = 1e-3)    # Poisson (Borel) offspring
  expect_equal(.chain_logpmf(3, 0.9, 50), .chain_logpmf(3, 0.9, Inf), tolerance = 0.02)  # NB -> Borel
})

test_that("cluster inference jointly recovers R and the dispersion k from chain sizes", {
  # a superspreading cluster set (k = 0.4) with enough introductions to pin both parameters
  cfg <- default_config(); cfg$clusters$dispersion_k <- 0.4; cfg$clusters$n_introductions <- 300
  set.seed(3); cl <- simulate_clusters(cfg)
  f <- cluster_size_fit(cl$true_sizes)
  expect_equal(f$R, 0.9, tolerance = 0.12)             # true R = 0.9
  expect_equal(f$k, 0.4, tolerance = 0.25)             # true dispersion k = 0.4
  expect_true(0.9 >= f$R_ci[1] && 0.9 <= f$R_ci[2])    # truth inside the CI
  expect_true(0.4 >= f$k_ci[1] && 0.4 <= f$k_ci[2])
})

test_that("the cluster tool runs on the observed sizes and scores against truth", {
  cf <- cluster_analysis(test_input)
  sc <- score_clusters(test_sim, cf)
  expect_gt(cf$R, 0)                                    # a positive reproduction number
  expect_gt(cf$k, 0)                                    # a positive dispersion
  expect_true(sc$R_in_ci)                               # the true R lies in the CI
  expect_true(is.finite(cf$extinction_prob))
})

test_that("validate_config rejects bad cluster parameters", {
  cfg <- default_config()
  bad <- cfg; bad$clusters$dispersion_k <- 0
  expect_error(validate_config(bad), "dispersion_k")
  bad2 <- cfg; bad2$clusters$detection_prob <- 1.5
  expect_error(validate_config(bad2), "detection_prob")
})

test_that("catchment back-calculation recovers the true source prevalence", {
  cb <- catchment_backcalc(test_input, 30:60, min_surveillance = 0.7,
                           source_pop = test_sim$config$source$population)
  sc <- score_catchment(test_sim, cb, test_input)
  expect_equal(sc$ratio, 1, tolerance = 0.2)                        # within 20% of truth
  expect_gt(sc$true_underascertainment, 2)                          # source really is much bigger than reported
})

test_that("visitor-vs-resident factor is a growing fraction in (0, 1]", {
  expect_equal(visitor_resident_factor(0, 0.15, 0.2), 0)               # zero stay -> no exposure window
  expect_lt(visitor_resident_factor(3, 0.15, 0.2),
            visitor_resident_factor(30, 0.15, 0.2))                    # longer stays are more representative
  expect_lt(visitor_resident_factor(30, 0.15, 0.2), 1)                 # never quite reaches 1
})

test_that("catchment range brackets the uncorrected estimate and is correctly ordered", {
  rg <- catchment_range(test_input, 30:60, source_pop = test_sim$config$source$population,
                        growth_rate = 0.15, recovery_rate = 0.2, min_surveillance = 0.7)
  cb <- catchment_backcalc(test_input, 30:60, min_surveillance = 0.7,
                           source_pop = test_sim$config$source$population)
  expect_equal(rg$low, cb$est_prevalence, tolerance = 1e-6)            # low end == uncorrected floor
  expect_lte(rg$low, rg$central)                                       # corrections only ever push UP
  expect_lte(rg$central, rg$high)
})

test_that("importation-risk flagging tracks true surveillance quality", {
  ir <- importation_risk(test_input, window = 25:70, min_surveillance = 0.7)
  sc <- score_importation_risk(ir)
  expect_gt(sc$cor_residual_surveillance, 0.4)                      # residuals rise with true quality
  if (!is.na(sc$mean_surv_flagged) && !is.na(sc$mean_surv_unflagged))
    expect_lt(sc$mean_surv_flagged, sc$mean_surv_unflagged)         # flagged countries are the poorly surveilled ones
})

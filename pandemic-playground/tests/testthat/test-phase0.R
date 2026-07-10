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

test_that("catchment back-calculation recovers the true source prevalence", {
  cb <- catchment_backcalc(test_input, 30:60, min_surveillance = 0.7,
                           source_pop = test_sim$config$source$population)
  sc <- score_catchment(test_sim, cb, test_input)
  expect_equal(sc$ratio, 1, tolerance = 0.2)                        # within 20% of truth
  expect_gt(sc$true_underascertainment, 2)                          # source really is much bigger than reported
})

test_that("importation-risk flagging tracks true surveillance quality", {
  ir <- importation_risk(test_input, window = 25:70, min_surveillance = 0.7)
  sc <- score_importation_risk(ir)
  expect_gt(sc$cor_residual_surveillance, 0.4)                      # residuals rise with true quality
  if (!is.na(sc$mean_surv_flagged) && !is.na(sc$mean_surv_unflagged))
    expect_lt(sc$mean_surv_flagged, sc$mean_surv_unflagged)         # flagged countries are the poorly surveilled ones
})

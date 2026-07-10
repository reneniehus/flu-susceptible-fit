# Phase 3 tools (analysis/phase3_*.R).

test_that("the SIRS integrator reproduces the analytic SIR final-size relation", {
  for (R0 in c(1.3, 1.8, 2.5)) {
    tr <- sirs_integrate(R0 = R0, gamma = 1 / 5, omega = 0, S0 = 1 - 1e-4, I0 = 1e-4, days = 600, dt = 0.1)
    z_sim    <- 1 - tail(tr$S, 1)
    z_theory <- uniroot(function(z) 1 - exp(-R0 * z) - z, c(1e-6, 1 - 1e-9))$root
    expect_equal(z_sim, z_theory, tolerance = 1e-3)
  }
})

test_that("SIRS scenarios: a fitter variant raises the peak, and boosting lowers it", {
  grid <- list(waning = c(slow = 1/365, fast = 1/120),
               R0 = c(wildtype = 1.3, variant = 1.8), uptake = c(none = 0, campaign = 0.4))
  res <- sirs_scenarios(grid, immune0 = 0.1, days = 365)
  variant_peak  <- mean(res$peak_prevalence[res$R0 == "variant"  & res$uptake == "none"])
  wildtype_peak <- mean(res$peak_prevalence[res$R0 == "wildtype" & res$uptake == "none"])
  expect_gt(variant_peak, wildtype_peak)                            # variant is worse
  be <- boosting_effect(res, "variant")
  expect_gt(be$peak_reduction, 0)                                   # boosting helps
})

test_that("the selection coefficient recovers the realized variant advantage", {
  vs <- variant_selection(test_input, "IT", window = 40:150)
  sc <- score_variant_selection(test_sim, vs, window = 40:150)
  expect_gt(vs$s, 0)                                                # variant is advancing
  expect_equal(vs$s, sc$s_realized, tolerance = 0.03)              # matches the realized logit slope
})

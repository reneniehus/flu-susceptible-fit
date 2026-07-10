# Phase 2 tools (analysis/phase2_*.R).

test_that("renewal_project reproduces a known exponential growth rate", {
  gi <- discretise(epidist_gamma("gi", 5, 1.7), boundary = "cori")
  hist <- exp(0.1 * (0:40))                                         # a clean exponential seed
  proj <- renewal_project(hist, R = r_to_R(0.1, gi), gi, horizon = 20)
  r_hat <- mean(diff(log(proj)))
  expect_equal(r_hat, 0.1, tolerance = 0.02)                        # projects at ~the intended rate
})

test_that("the capacity forecast calls a breach correctly in the clear growth phase", {
  cap <- truth_capacity(test_sim, "IT")
  fc <- forecast_capacity(test_input, "IT", as_of = 90, horizon = 28, capacity = cap)
  sc <- score_forecast(test_sim, fc, test_input)
  expect_true(sc$forecast_breach)                                   # forecast sees the breach coming
  expect_true(sc$true_breach)                                       # and the truth breaches
  expect_lt(abs(sc$forecast_breach_day - sc$true_breach_day), 10)   # timing within ~10 days
})

test_that("the intervention ITS detects the imposed slowing of transmission", {
  for (cc in c("IT", "PL", "DE")) {
    iday <- test_sim$config$rt_country[[cc]]$t[2]
    it <- intervention_its(test_input, cc, intervention_day = iday, window = 21)
    sc <- score_intervention(test_sim, it, window = 21)
    expect_true(it$slowed)                                          # significant slowing detected
    expect_lt(it$R_after, it$R_before)                             # R drops
    expect_true(sc$true_slowing)                                    # and it really did slow
  }
})

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

test_that("the capacity forecast stays physical (bounded by population) for a fast, high-R pathogen", {
  # a fast pathogen where a naive fixed-R projection (no depletion) would explode past the population
  cfg <- config_subset(default_config(), c("IT", "DE", "PL"), n_days = 120)
  cfg$delays$generation_interval <- epidist_gamma("gi", 2.6, 1.3)
  cfg$rt_source  <- step_schedule(c(0, 55), c(2.0, 0.9))
  cfg$rt_country <- make_country_rt(cfg$countries, R_start = 2.0, R_post = 0.8, intervention_days = c(60, 65, 70))
  sim <- simulate_pandemic(cfg); inp <- as_analysis_input(sim)
  pop <- sim$par$country_pop[["IT"]]
  fc  <- forecast_capacity(inp, "IT", as_of = 50, horizon = 28,
                           capacity = truth_capacity(sim, "IT"), population = pop)
  capped   <- max(fc$projection$admissions, na.rm = TRUE)
  uncapped <- max(forecast_capacity(inp, "IT", as_of = 50, horizon = 28,
                                    capacity = truth_capacity(sim, "IT"), population = NULL)$projection$admissions,
                  na.rm = TRUE)
  expect_true(is.finite(capped))
  expect_lt(capped, pop)                               # cannot imply more admissions than the population
  expect_lt(capped, uncapped)                          # depletion/cap really bit (the naive forecast is larger)
})

test_that("the forecast flags an out-of-envelope projection that shoots through capacity", {
  fc   <- forecast_capacity(test_input, "IT", as_of = 90, horizon = 28, capacity = truth_capacity(test_sim, "IT"))
  peak <- max(fc$breach$peak_admissions)
  # a capacity far BELOW the projected peak: the curve blasts through it -> out of envelope (unphysical)
  fc_low  <- forecast_capacity(test_input, "IT", as_of = 90, horizon = 28, capacity = peak / 10)
  expect_true(fc_low$out_of_envelope)
  expect_true(any(fc_low$breach$out_of_envelope))
  # a capacity comfortably ABOVE the peak: never approached -> in envelope
  fc_high <- forecast_capacity(test_input, "IT", as_of = 90, horizon = 28, capacity = peak * 2)
  expect_false(fc_high$out_of_envelope)
})

test_that("intervention_its exposes a depletion-suspected flag", {
  iday <- test_sim$config$rt_country[["IT"]]$t[2]
  it <- intervention_its(test_input, "IT", intervention_day = iday, window = 21)
  expect_type(it$depletion_suspected, "logical")
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

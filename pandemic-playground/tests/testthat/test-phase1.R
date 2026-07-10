# Phase 1 tools recover the truth (analysis/phase1_*.R).

test_that("delay-adjusted CFR exceeds naive CFR early and tracks the true confirmed CFR", {
  s70  <- cfr_static(test_input, "IT", as_of = 70)
  s150 <- cfr_static(test_input, "IT", as_of = 150)
  expect_gt(s70$cfr_adjusted, s70$cfr_naive)                        # adjustment lifts the early estimate
  sk <- score_cfr(test_sim, s150, test_input)
  expect_equal(s150$cfr_adjusted, sk$true_confirmed_cfr, tolerance = 0.4 * sk$true_confirmed_cfr)
  expect_gt(sk$ratio_cfr_to_ifr, 1)                                 # CFR overstates IFR (case under-ascertainment)
})

test_that("Cori Rt recovers the realized Rt on the (import-free) source infection curve", {
  # the source X has no importation, so the renewal assumption holds exactly and Cori recovers it well
  gi <- discretise(test_input$delays$generation_interval, boundary = "cori")
  inf <- truth_infections(test_sim, "X")
  est <- estimate_rt_cori(inf$infections, inf$day, gi, window = 7); est$location <- "X"
  sc <- score_rt(test_sim, est)
  expect_gt(sc$cor, 0.9)                                            # tracks the true trajectory closely
  expect_lt(sc$mae, 0.2)
})

test_that("importation inflates a country's EARLY Cori Rt above the endogenous realized Rt", {
  # a documented limitation: Cori attributes imported infections to local transmission, so early on
  # (when imports are a large share of local incidence) it reads high. This is a feature to be aware of.
  gi <- discretise(test_input$delays$generation_interval, boundary = "cori")
  inf <- truth_infections(test_sim, "IT")
  est <- estimate_rt_cori(inf$infections, inf$day, gi, window = 7)
  rt  <- truth_rt(test_sim, "IT")
  early <- est$day %in% 25:45
  est_early  <- mean(est$Rt[early])
  true_early <- mean(rt$Rt_effective[match(est$day[early], rt$day)])
  expect_gt(est_early, true_early)                                 # biased high while importation dominates
})

test_that("Cori Rt on true infections correctly crosses 1 at the intervention", {
  gi <- discretise(test_input$delays$generation_interval, boundary = "cori")
  inf <- truth_infections(test_sim, "IT")
  est <- estimate_rt_cori(inf$infections, inf$day, gi, window = 7)
  iday <- test_sim$config$rt_country[["IT"]]$t[2]
  before <- est$Rt[est$day %in% (iday - 20):(iday - 5)]
  after  <- est$Rt[est$day %in% (iday + 10):(iday + 25)]
  expect_gt(mean(before), 1)                                        # growing before
  expect_lt(mean(after), 1)                                         # controlled after
})

test_that("the nowcast beats the raw truncated counts on recent onset days", {
  nc <- nowcast_truncation(test_input, "IT", as_of = 110)
  sc <- score_nowcast(test_sim, nc, recent = 14)
  expect_gt(sc$rmse_improvement, 0.3)                              # at least 30% better than not nowcasting
  # the leading-edge honesty guard flags the most recent (near-zero-completeness) day(s)
  expect_true(any(nc$flagged_too_recent))
  # the improvement is a LIKE-FOR-LIKE comparison: both scored on the same day set
  expect_equal(sc$nowcast$n, sc$observed_naive$n)
})

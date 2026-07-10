# The observation model degrades the truth in the right, quantifiable ways (R/observe.R).

test_that("deaths are approximately IFR x death-detection of infections", {
  # over the whole run, total detected deaths ~ IFR * death_rho * total infections (delays only shift timing)
  inf   <- sum(loc_series(test_sim$truth$infections, "IT")$infections)
  death <- sum(loc_series(test_sim$observed$deaths, "IT")$deaths_by_date)
  expected <- inf * test_sim$par$ifr_eff * test_sim$config$ascertainment$death_rho
  # allow right-censoring of the latest deaths + Poisson noise
  expect_gt(death, 0.7 * expected)
  expect_lt(death, 1.05 * expected)
})

test_that("case ascertainment thins onsets, and rises over time", {
  # ratio of reported cases to true onsets should be well below 1 early (rho ~ 0.05) and higher later
  loc <- "IT"
  onsets <- loc_series(test_sim$truth$onsets, loc)
  cases  <- loc_series(test_sim$truth$cases_by_onset, loc)          # eventual detected (not truncated)
  early <- onsets$day %in% 20:35; late <- onsets$day %in% 80:100
  r_early <- sum(cases$cases[early]) / max(sum(onsets$onsets[early]), 1)
  r_late  <- sum(cases$cases[late])  / max(sum(onsets$onsets[late]),  1)
  expect_lt(r_early, r_late)                                         # ascertainment ramps up
  expect_lt(r_early, 0.3)                                            # early detection is poor
})

test_that("detected imports are a thinned copy of true imports (never more)", {
  imp_true <- test_sim$truth$imports
  imp_obs  <- test_sim$observed$detected_imports
  for (cc in test_sim$config$countries) {
    t <- sum(imp_true$imports[imp_true$country == cc])
    o <- sum(imp_obs$detected_imports[imp_obs$country == cc])
    expect_lte(o, t)                                                # detection cannot exceed truth
  }
  # a high-surveillance country detects a larger FRACTION than a low-surveillance one
  sq <- test_sim$config$surveillance_quality
  hi <- names(which.max(sq)); lo <- names(which.min(sq))
  frac <- function(cc) sum(imp_obs$detected_imports[imp_obs$country == cc]) /
                       max(sum(imp_true$imports[imp_true$country == cc]), 1)
  expect_gt(frac(hi), frac(lo))
})

test_that("observed flight volumes are a noisy, unbiased copy of the true volumes", {
  vt <- test_sim$truth$flight_volumes; vo <- test_sim$observed$flight_volumes
  m <- merge(vt, vo, by = c("country", "day"), suffixes = c("_true", "_obs"))
  expect_gt(stats::cor(m$volume_true, m$volume_obs), 0.98)          # tracks closely
  expect_equal(mean(m$volume_obs) / mean(m$volume_true), 1, tolerance = 0.03)   # ~unbiased
})

test_that("admissions accumulate and are dated after infections", {
  loc <- "IT"
  adm <- loc_series(test_sim$observed$admissions, loc)
  inf <- loc_series(test_sim$truth$infections, loc)
  expect_gt(sum(adm$admissions), 0)
  expect_gt(which.max(adm$admissions), which.max(inf$infections))    # admissions peak after infections
})

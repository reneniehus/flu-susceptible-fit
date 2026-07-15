# The full pipeline: reproducibility, structure, and the truth/observed contract (R/assemble.R, observe.R).

test_that("simulate_pandemic returns the two-worlds structure", {
  expect_s3_class(test_sim, "pandemic_sim")
  expect_true(all(c("truth", "observed", "config", "par", "latent") %in% names(test_sim)))
  expect_true(all(c("infections", "Rt", "cases_by_onset", "imports", "flight_volumes",
                    "variant_freq", "capacity", "prevalence") %in% names(test_sim$truth)))
  expect_true(all(c("cases_by_onset", "cases_by_report", "reporting_triangle", "deaths",
                    "admissions", "detected_imports", "variant_cases", "flight_volumes") %in%
                    names(test_sim$observed)))
})

test_that("the run is reproducible: same config + seed -> identical data", {
  a <- simulate_pandemic(test_cfg)
  b <- simulate_pandemic(test_cfg)
  expect_identical(a$observed$cases_by_report, b$observed$cases_by_report)
  expect_identical(a$observed$deaths, b$observed$deaths)
  expect_identical(a$truth$infections, b$truth$infections)
  # a different seed gives different data
  cfg2 <- test_cfg; cfg2$seed <- test_cfg$seed + 1
  expect_false(isTRUE(all.equal(simulate_pandemic(cfg2)$observed$deaths, a$observed$deaths)))
})

test_that("default_config() is deterministic yet leaves the caller's RNG untouched", {
  a <- default_config(); b <- default_config()
  expect_identical(a, b)                                # the default config is byte-identical every call
  set.seed(42); x <- stats::runif(1)
  set.seed(42); invisible(default_config()); y <- stats::runif(1)
  expect_identical(x, y)                                # building a config must not disturb the RNG stream
})

test_that("the source epidemic seeds the countries by importation", {
  imp <- test_sim$truth$imports
  seeded <- tapply(imp$imports, imp$country, sum)
  expect_true(all(seeded > 0))                                       # every country receives imports
  # bigger, better-connected countries import more (DE >> DK here)
  expect_gt(seeded[["DE"]], seeded[["DK"]])
})

test_that("observed cases-by-onset are right-truncated relative to the eventual truth", {
  for (loc in c("X", "IT")) {
    obs <- loc_series(test_sim$observed$cases_by_onset, loc)$cases
    evt <- loc_series(test_sim$truth$cases_by_onset, loc)$cases
    expect_true(all(obs <= evt + 1e-9))                              # observed never exceeds eventual
    tail_obs <- tail(obs, 10); tail_evt <- tail(evt, 10)
    expect_lt(sum(tail_obs), sum(tail_evt))                          # recent onsets visibly undercounted
  }
})

test_that("the reporting triangle reconciles with the by-report curve", {
  loc <- "IT"
  tri <- test_sim$observed$reporting_triangle
  rep_from_tri <- tapply(tri$cases[tri$location == loc], tri$report_day[tri$location == loc], sum)
  rep_curve    <- loc_series(test_sim$observed$cases_by_report, loc)
  total_tri <- sum(rep_from_tri); total_curve <- sum(rep_curve$cases)
  expect_equal(total_tri, total_curve)                              # every within-horizon report is counted once
})

test_that("turning the variant off yields a single-strain run with zero variant frequency", {
  cfg <- test_cfg; cfg$variant$enabled <- FALSE
  sim <- simulate_pandemic(cfg)
  vf <- loc_series(sim$truth$variant_freq, "X")$variant_freq
  expect_true(all(vf == 0))
  expect_equal(sum(sim$observed$variant_cases$variant), 0)
})

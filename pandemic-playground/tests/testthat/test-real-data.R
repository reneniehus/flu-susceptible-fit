# The real-data seam: tools consume the tidy SCHEMA, not the sim object, so hand-built ("real") data
# in the same shape works identically (analysis/analysis_common.R).

test_that("as_analysis_input exposes the documented schema", {
  expect_true(all(c("cases_by_onset", "cases_by_report", "reporting_triangle", "deaths", "admissions",
                    "detected_imports", "variant_cases", "flight_volumes", "delays", "as_of",
                    "countries", "source_code") %in% names(test_input)))
  expect_true(inherits(test_input$delays$generation_interval, "epidist"))
})

test_that("a tool gives the same answer whether fed the sim adapter or a hand-built schema list", {
  # a real-data user assembles this list by hand from their own surveillance tables
  hand <- list(
    cases_by_report = test_sim$observed$cases_by_report,
    cases_by_onset  = test_sim$observed$cases_by_onset,
    delays          = list(generation_interval = epidist_gamma("gi", 5.2, 1.7)),
    source_code     = "X"
  )
  a <- growth_analysis(test_input, "X", window = 20:45)
  b <- growth_analysis(hand,       "X", window = 20:45)
  expect_equal(a$r, b$r)
  expect_equal(a$R, b$R)
})

test_that("the growth tool recovers a known rate from a purely hand-made dataset", {
  # no simulator at all: a clean exponential case series a user might supply
  day <- 0:40; r_true <- 0.12
  cases <- rpois(length(day), 5 * exp(r_true * day))
  real <- list(cases_by_report = data.frame(location = "REAL", day = day, cases = cases),
               delays = list(generation_interval = epidist_gamma("gi", 5, 1.7)),
               source_code = "REAL")
  ga <- growth_analysis(real, "REAL", window = day)
  expect_equal(ga$r, r_true, tolerance = 0.03)                      # recovers the planted growth rate
})

test_that("variant selection runs on a hand-built sequencing table", {
  day <- seq(0, 120, by = 3); s_true <- 0.06
  freq <- plogis(-4 + s_true * day)
  seqd <- rep(200L, length(day)); variant <- rbinom(length(day), seqd, freq)
  real <- list(variant_cases = data.frame(location = "REAL", day = day, sequenced = seqd, variant = variant))
  vs <- variant_selection(real, "REAL")
  expect_equal(vs$s, s_true, tolerance = 0.02)                      # recovers the planted selection coefficient
})

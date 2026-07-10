# Delay distributions and their discretisation (R/epidist.R).

test_that("convenience constructors recover the requested mean", {
  expect_equal(epidist_mean(epidist_gamma("g", 5.2, 1.7)), 5.2, tolerance = 1e-8)
  expect_equal(epidist_mean(epidist_lognormal("ln", 5.5, 2.3)), 5.5, tolerance = 1e-8)
})

test_that("discretisation is a proper PMF whose mean matches the distribution", {
  d <- epidist_gamma("g", 5.2, 1.7); w <- discretise(d)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(all(w >= 0))
  expect_equal(sum((0:(length(w) - 1)) * w), 5.2, tolerance = 0.05)   # small discretisation error
})

test_that("the 'cori' boundary zeros the day-0 (no same-day) generation-interval mass", {
  gi <- epidist_gamma("gi", 5.2, 1.7)
  expect_equal(discretise(gi, boundary = "cori")[1], 0)
  expect_gt(discretise(gi, boundary = "interval")[1], 0)             # incubation-style keeps day 0
})

test_that("convolving two delay PMFs adds their means (used for infection->death)", {
  a <- discretise(epidist_lognormal("inc", 5.5, 2.3))
  b <- discretise(epidist_gamma("o2d", 15, 6.6))
  conv <- stats::convolve(a, rev(b), type = "open"); conv[conv < 0] <- 0; conv <- conv / sum(conv)
  m <- function(p) sum((0:(length(p) - 1)) * p)
  expect_equal(m(conv), m(a) + m(b), tolerance = 0.1)
})

test_that("epidist construction validates its parameters", {
  expect_error(epidist("bad", "gamma", list(shape = 2)))            # missing rate
  expect_error(epidist("bad", "nonsense", list()))                 # unknown family
})

test_that("convenience constructors reject non-positive mean/sd", {
  expect_error(epidist_gamma("x", mean = 5, sd = 0))               # sd = 0 -> Inf params
  expect_error(epidist_gamma("x", mean = -5, sd = 2))              # negative mean
  expect_error(epidist_lognormal("x", mean = 5, sd = 0))
  expect_error(epidist_weibull("x", shape = 0, scale = 2))
})

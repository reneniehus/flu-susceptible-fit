# utils.R
#
# Small, dependency-light helpers shared across the simulator and the analysis toolbox: evaluating
# step-function schedules, laying out the day/date grid, seasonality, and the causal delay
# convolution that turns "events on day t" into "downstream events on later days" (infections->onsets,
# onsets->deaths, onsets->reports). Base R only; each helper does one thing.

# ---- |-evaluate a change-point (step) schedule at day(s) t ----
# schedule = list(t = change-point days, value = level from each change-point). Vectorised over t.
# `step_at(list(t=c(0,40,70), value=c(2.4,1.2,0.9)), c(10, 50, 80))` -> c(2.4, 1.2, 0.9).
step_at <- function(schedule, t) {
  idx <- findInterval(t, schedule$t)                 # which segment each t falls in (t[1] must be 0)
  schedule$value[pmax(idx, 1)]
}

# ---- |-the day / date grid for a run (day index is 0-based) ----
# Everything internal is indexed by integer day 0..(n_days-1); `date` carries the calendar mapping
# for outputs and for turning config dates (travel-ban, variant intro) into day indices.
sim_days <- function(cfg) {
  day  <- 0:(cfg$n_days - 1)
  list(day = day, date = cfg$start_date + day, n = cfg$n_days)
}

# ---- |-turn a calendar date into a 0-based day index within the run (NA if outside) ----
date_to_day <- function(date, cfg) {
  d <- as.integer(as.Date(date) - cfg$start_date)
  if (is.na(d) || d < 0 || d >= cfg$n_days) return(NA_integer_)
  d
}

# ---- |-mild annual seasonality multiplier (>= 0), peaking on `peak_day` ----
# Used to modulate flight volumes over the year: 1 + amp*cos(2pi (day - peak_day)/period).
seasonality_multiplier <- function(day, amp, peak_day, period = 365) {
  pmax(0, 1 + amp * cos(2 * pi * (day - peak_day) / period))
}

# ---- |-causal delay convolution: place day-t events onto later days via a delay PMF ----
# Given a daily series x (x[t] events on day t, t = 1..n for days 0..n-1) and a delay PMF `pmf`
# indexed from delay 0 (pmf[1] = P(delay 0)), returns y with y[t] = sum_{d>=0} x[t-d] * pmf[d+1].
# The tail that would fall beyond day n-1 is dropped -- this is exactly the right-censoring an analyst
# faces (deaths / reports not yet observed), and is what makes the observation model realistic.
shift_by_delay <- function(x, pmf) {
  n <- length(x); D <- length(pmf) - 1
  y <- numeric(n)
  for (d in 0:D) {
    if (pmf[d + 1] == 0) next
    # contribution of x to days (d+1)..n, i.e. shifted forward by d days
    if (d < n) y[(d + 1):n] <- y[(d + 1):n] + pmf[d + 1] * x[1:(n - d)]
  }
  y
}

# ---- |-full (uncensored) delay convolution, extending the output by the delay support ----
# Same as shift_by_delay() but keeps the tail beyond day n-1, returning a longer vector. Used for the
# LATENT truth (e.g. all deaths eventually caused), where we do not want to censor. Length n + D.
shift_by_delay_full <- function(x, pmf) {
  n <- length(x); D <- length(pmf) - 1
  y <- numeric(n + D)
  for (d in 0:D) {
    if (pmf[d + 1] == 0) next
    y[(d + 1):(n + d)] <- y[(d + 1):(n + d)] + pmf[d + 1] * x
  }
  y
}

# ---- |-lognormal multiplicative measurement noise with a given coefficient of variation ----
# For OBSERVED flight volumes: multiply true volumes by exp(N(-s^2/2, s^2)) so the noise is unbiased
# on the natural scale (median-unbiased, mean-preserving) with the requested CV.
rlnorm_measure <- function(x, cv) {
  s <- sqrt(log1p(cv^2))
  x * exp(stats::rnorm(length(x), -s^2 / 2, s))
}

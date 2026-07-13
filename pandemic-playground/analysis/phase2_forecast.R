# phase2_forecast.R
#
# PHASE 2 -- established transmission, growth to peak, healthcare demand.
# QUESTION: "Will we breach hospital / ICU capacity in the next few weeks?"
#
# THE LEAN METHOD (renewal forecast of cases -> admissions, holding R fixed). Admissions are a delayed,
# scaled copy of infections, so we forecast the thing that LEADS them -- the case curve -- and then map
# it through the onset->admission delay:
#   1. reconstruct the recent case curve, nowcast-correcting its right-truncated tail (Phase 1);
#   2. estimate the current reproduction number R from those cases (Cori);
#   3. project cases forward by a fixed-R renewal, over a small SCENARIO SET (e.g. R, 0.8R, 1.2R);
#   4. map the full (observed + projected) case trajectory to admissions via a fitted case->admission
#      ratio and the onset->admission delay;  5. compare to capacity to read off IF and WHEN it breaks.
#
# WHY CASES, NOT ADMISSIONS, DRIVE THE FORECAST. Estimating R off admissions is lagged by the
# admission delay, so near the peak it still reads "growing" when transmission has already turned over,
# and the forecast overshoots. Cases lead admissions, so the case-driven forecast peaks on time and the
# admission delay is put back in explicitly -- admissions keep rising for a week or two after cases
# peak (the demand already "baked in"), exactly as in reality. Holding R fixed still cannot predict a
# NEW change in transmission -- that is what the scenario set is for, and it is stated as association.
#
# THE ONE LIMIT TO STATE PLAINLY. R here is estimated from ONSET-dated cases, so it lags the true
# INFECTION-time R by about one incubation period: when transmission drops sharply (a lockdown), the
# forecast keeps reading the old, higher R for several days and OVERSHOOTS right at the turnover. The
# playground shows this by scoring a forecast made at the peak against what actually happened.
# Deconvolving to infection dates first (EpiNow2) removes most of the lag -- a documented extension.
#
# References
#   Held L, Meyer S, Bracher J. Probabilistic forecasting in infectious disease epidemiology. Stat Med. 2017.
#   EpiNow2 (Epiverse-TRACE) -- the fuller renewal forecast (infers infections; uncertainty end-to-end).
#
# Requires: analysis_common.R, phase1_rt.R, phase1_nowcast.R; R/epidist.R (discretise).

# ---- |-project a series forward by a fixed-R renewal step, with susceptible depletion ----
# future_t = R * (S_t / pool) * sum_s gi[s] * series_{t-s}, extending `history` day by day. `pool` is
# the remaining susceptibles on the SERIES scale (Inf = no depletion). Without depletion a fixed R > 1
# projects UNBOUNDED exponential growth, so for a fast / high-R pathogen the naive forecast implied
# more admissions than the entire population; the S/pool factor turns the projection over, and the cap
# guarantees cumulative projected incidence never exceeds the susceptible pool. Returns just the tail.
renewal_project <- function(history, R, gi_pmf, horizon, pool = Inf) {
  D <- length(gi_pmf) - 1; series <- as.numeric(history); n0 <- length(series)
  remaining <- pool
  for (h in seq_len(horizon)) {
    t <- n0 + h; smax <- min(D, t - 1)
    force <- sum(gi_pmf[2:(smax + 1)] * series[(t - 1):(t - smax)])
    depletion <- if (is.finite(pool)) max(remaining / pool, 0) else 1     # remaining susceptible fraction
    projected <- R * depletion * force
    if (is.finite(pool)) projected <- min(projected, max(remaining, 0))   # cannot exceed the pool
    series[t] <- projected; remaining <- remaining - projected
  }
  series[(n0 + 1):(n0 + horizon)]
}

# ---- |-forecast admissions from the case curve and test them against capacity, over R scenarios ----
# input : as_analysis_input(sim). location, as_of, horizon, capacity as before. rt : optional fixed R
# (else estimated from recent cases). scenarios : multiplicative R factors. nowcast : correct the
# truncated case tail before forecasting. population + susceptible_fraction : the analyst's country
# population and (e.g. serology-based) remaining-susceptible fraction, used to deplete/cap the
# projection so it stays physical -- essential for fast/high-R pathogens (NULL population = no cap,
# fine only well below depletion). rt_window : NULL derives it from the GI (as in the Cori tool).
forecast_capacity <- function(input, location, as_of, horizon = 28, capacity,
                              rt = NULL, scenarios = c(0.8, 1.0, 1.2),
                              rt_window = NULL, nowcast = TRUE, gi = NULL,
                              population = NULL, susceptible_fraction = 1) {
  gi_dist <- gi %||% input$delays$generation_interval
  gi_pmf  <- discretise(gi_dist, boundary = "cori")
  o2a     <- discretise(input$delays$onset_to_admission)
  if (is.null(rt_window)) rt_window <- max(3L, round(1.5 * epidist_mean(gi_dist)))   # GI-scaled

  # ---- (1) the recent case curve, tail nowcast-corrected ----
  ca  <- loc_series(input$cases_by_onset, location); ca <- ca[ca$day <= as_of, ]
  cases <- ca$cases
  if (nowcast) {
    nc  <- nowcast_truncation(input, location, as_of = as_of)
    est <- nc$nowcast[match(ca$day, nc$onset_day)]        # nowcast estimate aligned to the case days
    use <- is.finite(est)                                 # use it where the nowcast produced a value
    cases[use] <- est[use]                                # else keep the raw (truncated) count
  }

  # ---- (2) current R from the (corrected) case curve ----
  R_now <- if (!is.null(rt)) rt else {
    est <- estimate_rt_cori(cases, ca$day, gi_pmf, window = rt_window)
    if (nrow(est)) stats::median(utils::tail(est$Rt, 5)) else NA_real_
  }

  # ---- fit the case->admission ratio on FULLY-OBSERVED days that actually contain admissions ----
  # A fixed "last 60 days" lookback contains zero admissions early in a fast epidemic -> chr = 0 -> a
  # silent zero forecast while a breach is underway. So fit on whatever settled days have admissions.
  fit_lookback <- 60L                                    # ~2 months of history to average the ratio over
  settle_lag   <- which(cumsum(o2a) >= 0.9)[1] - 1       # days for ~90% of a day's admissions to accrue
  adm      <- loc_series(input$admissions, location); adm <- adm$admissions[match(ca$day, adm$day)]
  case_adm <- shift_by_delay(cases, o2a)                 # cases mapped through the admission delay
  settled  <- ca$day <= (as_of - settle_lag) & ca$day >= (as_of - fit_lookback) & !is.na(adm) & adm > 0
  chr <- if (any(settled)) sum(adm[settled]) / max(sum(case_adm[settled]), 1e-9) else NA_real_

  # susceptible pool on the case scale: a valid (conservative) cap is the susceptible headcount, since
  # cases <= infections <= susceptibles. Supply serology via `susceptible_fraction` to sharpen it.
  pool <- if (!is.null(population)) susceptible_fraction * population else Inf

  # remaining susceptible pool = full pool minus what has already been infected (cases so far, a
  # conservative under-count on the case scale) -- so the projection starts from the right depletion.
  pool_remaining <- if (is.finite(pool)) max(pool - sum(cases), 0) else Inf

  # ---- (3-4) per scenario: project cases (with depletion), map to future admissions ----
  future_day <- (as_of + 1):(as_of + horizon)
  proj <- lapply(scenarios, function(s) {
    fut_cases  <- renewal_project(cases, R_now * s, gi_pmf, horizon, pool = pool_remaining)
    full_cases <- c(cases, fut_cases)                    # observed + projected
    full_adm   <- chr * shift_by_delay(full_cases, o2a)
    data.frame(day = future_day, scenario = sprintf("R x %.1f", s), R = R_now * s,
               admissions = full_adm[(length(cases) + 1):(length(cases) + horizon)])
  })
  proj <- do.call(rbind, proj)

  # ---- (5) breach day per scenario ----
  breach <- do.call(rbind, lapply(split(proj, proj$scenario), function(d) {
    hit <- which(d$admissions > capacity)
    data.frame(scenario = d$scenario[1], R = d$R[1], breaches = length(hit) > 0,
               breach_day = if (length(hit)) d$day[hit[1]] else NA_integer_,
               peak_admissions = max(d$admissions))
  }))

  list(location = location, as_of = as_of, horizon = horizon, R_now = R_now, chr = chr,
       capacity = capacity, projection = proj, breach = breach)
}

# ---- |-score the central forecast against the admissions that actually occurred ----
score_forecast <- function(sim, fc, input = as_analysis_input(sim)) {
  loc <- fc$location
  adm <- loc_series(input$admissions, loc)
  future <- fc$as_of + seq_len(fc$horizon)
  true_adm <- adm$admissions[match(future, adm$day)]

  central <- fc$projection[fc$projection$scenario == "R x 1.0", ]
  s <- pp_score(central$admissions, true_adm)

  hit <- which(true_adm > fc$capacity)                 # future days the truth actually exceeded capacity
  true_breach_day <- if (length(hit)) future[hit[1]] else NA_integer_
  cb <- fc$breach[fc$breach$scenario == "R x 1.0", ]
  list(location = loc, score = s,
       true_breach = !is.na(true_breach_day), true_breach_day = true_breach_day,
       forecast_breach = cb$breaches, forecast_breach_day = cb$breach_day)
}

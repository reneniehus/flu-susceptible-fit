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

# ---- |-project a series forward by a fixed-R renewal step (deterministic) ----
# future_t = R * sum_s gi[s] * series_{t-s}, extending `history` day by day. Returns just the tail.
renewal_project <- function(history, R, gi_pmf, horizon) {
  D <- length(gi_pmf) - 1; series <- as.numeric(history); n0 <- length(series)
  for (h in seq_len(horizon)) {
    t <- n0 + h; smax <- min(D, t - 1)
    series[t] <- R * sum(gi_pmf[2:(smax + 1)] * series[(t - 1):(t - smax)])
  }
  series[(n0 + 1):(n0 + horizon)]
}

# ---- |-forecast admissions from the case curve and test them against capacity, over R scenarios ----
# input : as_analysis_input(sim). location, as_of, horizon, capacity as before. rt : optional fixed R
# (else estimated from recent cases). scenarios : multiplicative R factors. nowcast : correct the
# truncated case tail before forecasting (recommended).
forecast_capacity <- function(input, location, as_of, horizon = 28, capacity,
                              rt = NULL, scenarios = c(0.8, 1.0, 1.2),
                              rt_window = 7, nowcast = TRUE, gi = NULL) {
  gi_pmf <- discretise(gi %||% input$delays$generation_interval, boundary = "cori")
  o2a    <- discretise(input$delays$onset_to_admission)

  # ---- (1) the recent case curve, tail nowcast-corrected ----
  ca  <- loc_series(input$cases_by_onset, location); ca <- ca[ca$day <= as_of, ]
  cases <- ca$cases
  if (nowcast) {
    nc <- nowcast_truncation(input, location, as_of = as_of)
    m  <- match(ca$day, nc$onset_day)
    fill <- ifelse(is.finite(nc$nowcast[m]), nc$nowcast[m], cases)
    cases <- ifelse(is.na(fill), cases, fill)
  }

  # ---- (2) current R from the (corrected) case curve ----
  R_now <- if (!is.null(rt)) rt else {
    est <- estimate_rt_cori(cases, ca$day, gi_pmf, window = rt_window)
    stats::median(utils::tail(est$Rt, 5))
  }

  # ---- fit the case->admission mapping (a scalar ratio) on a stable, fully-observed window ----
  adm <- loc_series(input$admissions, location); adm <- adm$admissions[match(ca$day, adm$day)]
  case_adm <- shift_by_delay(cases, o2a)                 # cases mapped through the admission delay
  fitwin <- ca$day <= (as_of - (length(o2a))) & ca$day >= (as_of - 60)   # avoid the truncated / delayed tail
  chr <- sum(adm[fitwin], na.rm = TRUE) / max(sum(case_adm[fitwin], na.rm = TRUE), 1e-9)

  # ---- (3-4) per scenario: project cases, map the full trajectory to future admissions ----
  future_day <- (as_of + 1):(as_of + horizon)
  proj <- lapply(scenarios, function(s) {
    fut_cases  <- renewal_project(cases, R_now * s, gi_pmf, horizon)
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

  true_breach_day <- { hit <- which(true_adm > fc$capacity); if (length(hit)) future[hit[1]] else NA_integer_ }
  cb <- fc$breach[fc$breach$scenario == "R x 1.0", ]
  list(location = loc, score = s,
       true_breach = !is.na(true_breach_day), true_breach_day = true_breach_day,
       forecast_breach = cb$breaches, forecast_breach_day = cb$breach_day)
}

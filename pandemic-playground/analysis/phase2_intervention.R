# phase2_intervention.R
#
# PHASE 2 -- established transmission.
# QUESTION: "Are our control measures working?"
#
# THE LEAN METHOD (interrupted time series on the growth rate). Fit a segmented log-linear (Poisson)
# model to the case curve with a break at the intervention date: the growth rate is r_before up to the
# break and r_after = r_before + delta afterwards. A negative delta -- growth slowing, the R crossing
# down through 1 -- is the signature of a working measure. We convert r_before / r_after to R through
# the generation interval (same Euler-Lotka map as Phase 0) so the effect reads as a before/after R.
#
# ASSOCIATION, NOT PROOF -- SAID PLAINLY. A step in the growth rate at the intervention date is
# consistent with the measure working, but it cannot PROVE it: seasonality, behaviour change,
# susceptible depletion and other simultaneous measures can all bend the curve at the same time. This
# tool quantifies an association and must be reported as one. (The playground can still SCORE it,
# because here we know the truth -- the R step really was imposed on that date.)
#
# THE LAG TO MIND. The break shows up in ONSET-dated cases about one incubation period AFTER the true
# intervention (and later still in report-dated cases). We therefore place the regression breakpoint at
# intervention_day + onset_lag so it lines up with where the effect actually appears in the data.
#
# References
#   Bernal JL, Cummins S, Gasparrini A. Interrupted time series regression for the evaluation of public
#     health interventions: a tutorial. Int J Epidemiol. 2017;46(1):348-355.
#
# Requires: analysis_common.R, phase0_growth_R.R (r_to_R); R/epidist.R (discretise).

# ---- |-segmented growth-rate (interrupted time series) around an intervention date ----
# input : as_analysis_input(sim). location, intervention_day. window : days each side to fit over.
# onset_lag : shift the breakpoint to where the effect appears in onset cases (default incubation mean).
intervention_its <- function(input, location, intervention_day, window = 21,
                             onset_lag = NULL, series = c("cases_by_onset", "cases_by_report"),
                             gi = NULL, level = 0.95) {
  series <- match.arg(series)
  gi_pmf <- discretise(gi %||% input$delays$generation_interval, boundary = "cori")
  if (is.null(onset_lag)) onset_lag <- round(epidist_mean(input$delays$incubation))
  bp <- intervention_day + onset_lag                    # the breakpoint as seen in the onset curve

  d <- loc_series(input[[series]], location)
  d <- d[d$day >= (bp - window) & d$day <= (bp + window), ]
  d <- d[is.finite(d$cases) & d$cases >= 0, ]
  if (nrow(d) < 8) stop("intervention_its: too few days around the breakpoint to fit")

  # hinge model: log E[cases] = a + r_before*(day-bp) + delta*max(day-bp, 0)
  t0 <- d$day - bp; hinge <- pmax(t0, 0)
  fit <- stats::glm(cases ~ t0 + hinge, family = stats::poisson(), data = d)
  co  <- stats::coef(fit); se <- sqrt(diag(stats::vcov(fit))); z <- stats::qnorm(1 - (1 - level) / 2)
  r_before <- unname(co["t0"]); r_after <- unname(co["t0"] + co["hinge"])
  delta    <- unname(co["hinge"]); delta_se <- unname(se["hinge"])

  list(location = location, intervention_day = intervention_day, breakpoint = bp,
       r_before = r_before, r_after = r_after,
       delta_r = delta, delta_ci = c(delta - z * delta_se, delta + z * delta_se),
       R_before = r_to_R(r_before, gi_pmf), R_after = r_to_R(r_after, gi_pmf),
       slowed = delta < 0 && (delta + z * delta_se) < 0)     # growth significantly slowed
}

# ---- |-score the estimated before/after R against the true realized R either side of the break ----
# Truth: realized R_t averaged over the fitting window before vs after the intervention. Confirms the
# tool detects the imposed step (and by how much it under/over-shoots it).
score_intervention <- function(sim, it, window = 21) {
  loc <- it$location; rt <- truth_rt(sim, loc)
  before <- rt$Rt_effective[rt$day >= (it$intervention_day - window) & rt$day < it$intervention_day]
  after  <- rt$Rt_effective[rt$day >  it$intervention_day & rt$day <= (it$intervention_day + window)]
  list(location = loc,
       R_before = it$R_before, true_R_before = mean(before, na.rm = TRUE),
       R_after  = it$R_after,  true_R_after  = mean(after,  na.rm = TRUE),
       detected_slowing = it$slowed,
       true_slowing = mean(after, na.rm = TRUE) < mean(before, na.rm = TRUE))
}

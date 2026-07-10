# phase1_nowcast.R
#
# PHASE 1 -- early local exponential growth.
# QUESTION: "How many cases really occurred recently, before the reports catch up?"
#
# THE LEAN METHOD. Standing at a data cutoff `as_of`, the cases already reported for a recent onset day
# o are only the fraction F(as_of - o) of the eventual total that report within the elapsed delay,
# where F is the reporting-delay CDF. Inverting that truncation gives the nowcast:
#     nowcast(o)  =  observed(o)  /  F(as_of - o)
# For onset days well in the past F -> 1 and nothing changes; for the last week or two F is small and
# the correction is large -- which is exactly where the raw onset curve misleads (it looks like the
# epidemic is turning over when reports simply have not arrived).
#
# TWO FLAVOURS. The delay CDF F can be the ASSUMED onset->report distribution (crude, needs no history)
# or ESTIMATED EMPIRICALLY from the triangle's completed onset days (defensible, self-correcting).
# Both are here. The fully-defensible version -- jointly modelling delays and reporting with
# uncertainty (epinowcast / EpiNow2::estimate_truncation) -- is a documented extension front.
#
# References
#   Hohle M, an der Heiden M. Bayesian nowcasting during the STEC O104:H4 outbreak. Biometrics. 2014.
#   epinowcast / EpiNow2 (Epiverse-TRACE) -- the fuller Bayesian nowcast this emulates.
#
# Requires: analysis_common.R; R/epidist.R (discretise) for the reporting delay.

# ---- |-empirical reporting-delay PMF from a triangle's completed onset days ----
# Uses only onset days old enough (as_of - onset >= max_delay) to be fully reported, so the estimated
# delay distribution is not itself truncated. Falls back to the assumed delay if too little history.
empirical_delay_pmf <- function(tri, as_of, max_delay) {
  tri <- tri[tri$report_day <= as_of, ]
  tri$delay <- tri$report_day - tri$onset_day
  complete <- tri[(as_of - tri$onset_day) >= max_delay, ]
  if (sum(complete$cases) < 50) return(NULL)              # not enough completed history -> caller uses assumed
  w <- tapply(complete$cases, factor(complete$delay, levels = 0:max_delay), sum)
  w[is.na(w)] <- 0
  as.numeric(w) / sum(w)
}

# ---- |-nowcast the recent onset curve by inverting the reporting-delay truncation ----
# input : as_analysis_input(sim). location : which location. as_of : the data cutoff (day index).
# delay_source : "assumed" (onset->report from config) or "empirical" (estimated from the triangle).
# Returns, for each onset day, the observed (truncated) count, the completeness F, and the nowcast
# with a predictive interval for the not-yet-reported tail.
nowcast_truncation <- function(input, location, as_of = input$as_of,
                               delay_source = c("assumed", "empirical"),
                               o2r = NULL, level = 0.95, min_completeness = 0.05) {
  delay_source <- match.arg(delay_source)
  tri <- input$reporting_triangle
  tri <- tri[tri$location == location & tri$report_day <= as_of, ]  # only what is visible at as_of

  assumed_pmf <- discretise(o2r %||% input$delays$onset_to_report)
  max_delay   <- length(assumed_pmf) - 1
  pmf <- if (delay_source == "empirical") empirical_delay_pmf(tri, as_of, max_delay) %||% assumed_pmf else assumed_pmf
  Fd  <- cumsum(pmf)

  # observed (truncated) cases by onset day, up to as_of
  onset_days <- 0:as_of
  observed <- tapply(tri$cases, factor(tri$onset_day, levels = onset_days), sum)
  observed[is.na(observed)] <- 0; observed <- as.numeric(observed)

  elapsed <- as_of - onset_days
  Fnow <- Fd[pmin(elapsed + 1, length(Fd))]              # completeness F(as_of - onset)
  point <- observed / Fnow

  # predictive interval for the not-yet-reported tail ~ Poisson(observed * (1-F)/F)
  lam <- observed * (1 - Fnow) / Fnow
  lo  <- (1 - level) / 2; hi <- 1 - lo
  lower <- observed + stats::qpois(lo, lam)
  upper <- observed + stats::qpois(hi, lam)

  # HONESTY GUARD: the leading edge (completeness below min_completeness) carries too little signal to
  # nowcast -- dividing by a near-zero F is unstable and the tail-Poisson interval is far too tight.
  # Flag those days NA rather than emit an overconfident number. This is a genuine limit of any
  # nowcast; the fuller Bayesian version narrows it but cannot remove it.
  too_recent <- Fnow < min_completeness
  point[too_recent] <- NA_real_; lower[too_recent] <- NA_real_; upper[too_recent] <- NA_real_

  data.frame(location = location, onset_day = onset_days, observed = observed,
             completeness = Fnow, nowcast = point, nowcast_lower = lower, nowcast_upper = upper,
             flagged_too_recent = too_recent)
}

# ---- |-score a nowcast against the eventual (fully-reported) onset totals ----
# Truth: the eventual detected cases by onset (sim$truth$cases_by_onset) -- what the truncated curve
# will converge to. We score the RECENT window (the last `recent` days before as_of), where truncation
# bites and the nowcast has to work; and compare it to just trusting the raw observed counts.
score_nowcast <- function(sim, nc, recent = 14) {
  loc <- nc$location[1]; as_of <- max(nc$onset_day)
  truth <- truth_cases_by_onset(sim, loc)
  Ntrue <- truth$cases[match(nc$onset_day, truth$day)]
  # LIKE-FOR-LIKE: score the nowcast and the naive observed on the SAME days -- the recent onset days
  # the nowcast actually produced an estimate for. The leading-edge days it flagged NA are excluded
  # from BOTH (neither has anything to say there), so the improvement is not inflated by crediting the
  # nowcast with the observed curve's error on days the nowcast declined to predict.
  win <- nc$onset_day > (as_of - recent) & is.finite(nc$nowcast)
  s_nowcast  <- pp_score(nc$nowcast[win],  Ntrue[win], nc$nowcast_lower[win], nc$nowcast_upper[win])
  s_observed <- pp_score(nc$observed[win], Ntrue[win])   # what you get by NOT nowcasting, SAME days
  list(location = loc, recent = recent, n_scored = s_nowcast$n,
       nowcast = s_nowcast, observed_naive = s_observed,
       rmse_improvement = 1 - s_nowcast$rmse / s_observed$rmse)
}

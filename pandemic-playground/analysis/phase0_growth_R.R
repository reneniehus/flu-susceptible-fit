# phase0_growth_R.R
#
# PHASE 0 -- before local introduction, at the source X.
# QUESTION: "Is it spreading person-to-person, and how fast?"
#
# THE LEAN METHOD. Fit exponential growth to the early case curve with a Poisson GLM (log-linear in
# time) to get the growth rate r; the doubling time is the back-of-envelope ln(2)/r. Then convert r to
# a reproduction number through the generation interval by the Euler-Lotka / Wallinga-Lipsitch
# identity: for a renewal process with generation-interval density g,
#     1 = R * integral g(a) e^{-r a} da        =>     R = 1 / sum_a g(a) e^{-r a}
# so R is r fed through the (discretised) generation interval. Growing epidemic <=> r > 0 <=> R > 1.
#
# WHY THIS AND NOT A MECHANISTIC FIT. Early on there is little data and strong right-truncation; a
# transparent log-linear growth rate plus a clearly-stated GI is more honest and more robust than a
# full renewal fit. The direct-R alternative (EpiEstim / a Cori renewal estimate) is implemented in
# phase1_rt.R and can be applied to the same early curve; here we keep the classic growth-rate route.
#
# References
#   Wallinga J, Lipsitch M. How generation intervals shape the relationship between growth rates and
#     reproductive numbers. Proc R Soc B. 2007;274(1609):599-604.
#   Wallinga J, Teunis P. Different epidemic curves for severe acute respiratory syndrome reveal
#     similar impacts of control measures. Am J Epidemiol. 2004;160(6):509-516. (related-work context)
#
# Requires: analysis_common.R; R/epidist.R (discretise) for the generation interval.

# ---- |-exponential growth rate from an early count series (Poisson GLM, log link) ----
# counts, day: aligned vectors over the fitting window. Returns the growth rate r (per day) with its
# 95% CI and the implied doubling time. A log-linear Poisson GLM is the standard early-growth fit.
estimate_growth_rate <- function(counts, day, level = 0.95) {
  ok <- is.finite(counts) & is.finite(day) & counts >= 0
  if (sum(ok) < 3) stop("estimate_growth_rate: need >= 3 points with data in the window")
  # too few events -> the Poisson GLM is unidentified and returns a meaningless r with an absurd CI;
  # refuse rather than emit a number that looks like an estimate (a barely-seeded country, say)
  if (sum(counts[ok]) < 5)
    stop("estimate_growth_rate: fewer than 5 events in the window -- too little signal to estimate growth")
  fit <- stats::glm(counts[ok] ~ day[ok], family = stats::poisson())
  r   <- unname(stats::coef(fit)[2])
  se  <- unname(sqrt(diag(stats::vcov(fit))[2]))
  if (!is.finite(se) || se > 5)                         # a near-degenerate fit: flag the huge CI honestly
    warning("estimate_growth_rate: the growth-rate estimate is poorly identified (very wide CI)")
  z   <- stats::qnorm(1 - (1 - level) / 2)
  list(r = r, r_lower = r - z * se, r_upper = r + z * se, se = se,
       doubling_time = log(2) / r,
       doubling_lower = log(2) / (r + z * se), doubling_upper = log(2) / (r - z * se))
}

# ---- |-convert a growth rate r to a reproduction number via the generation interval ----
# gi_pmf: discretised generation interval indexed from delay 0 (gi_pmf[1] = P(0) = 0). Vectorised
# over r, so a whole CI can be mapped through at once. This is the Euler-Lotka moment-generating
# transform R = 1 / M_g(-r) evaluated on the discrete GI.
r_to_R <- function(r, gi_pmf) {
  a <- 0:(length(gi_pmf) - 1)
  vapply(r, function(rr) 1 / sum(gi_pmf * exp(-rr * a)), numeric(1))
}

# ---- |-full Phase-0 growth analysis: r, doubling time and R (with CIs) from an early curve ----
# input : as_analysis_input(sim) (or a real-data list in the same schema).
# location : which location's case curve to fit (default the source X).
# window : day indices to fit over (the early, pre-turnover growth phase).
# series : "cases_by_report" (default) or "cases_by_onset"; gi : an epidist override for the GI.
growth_analysis <- function(input, location = input$source_code, window = 0:30,
                            series = c("cases_by_report", "cases_by_onset"), gi = NULL) {
  series <- match.arg(series)
  d      <- loc_series(input[[series]], location)
  d      <- d[d$day %in% window, ]
  gi_pmf <- discretise(gi %||% input$delays$generation_interval, boundary = "cori")

  g <- estimate_growth_rate(d$cases, d$day)
  R       <- r_to_R(g$r, gi_pmf)
  R_lower <- r_to_R(g$r_lower, gi_pmf)                # r_to_R is monotone increasing in r
  R_upper <- r_to_R(g$r_upper, gi_pmf)

  list(location = location, window = range(window), series = series,
       r = g$r, r_ci = c(g$r_lower, g$r_upper),
       doubling_time = g$doubling_time, doubling_ci = c(g$doubling_lower, g$doubling_upper),
       R = R, R_ci = c(R_lower, R_upper))
}

# ---- |-score a growth analysis against the true REALIZED R over the fitting window ----
# The Euler-Lotka r->R map recovers the R consistent with the observed growth rate -- i.e. the
# realized/effective R (strain-mix- and depletion-weighted), NOT the nominal schedule. So we score
# against `Rt_effective` by default (matching score_rt / score_intervention); the two coincide in the
# early single-strain window where growth analysis is normally used, but diverge once a window spans
# depletion or the variant. `which = "Rt"` scores against the nominal schedule instead.
score_growth <- function(sim, ga, which = c("Rt_effective", "Rt")) {
  which <- match.arg(which)
  rt   <- truth_rt(sim, ga$location)
  win  <- ga$window[1]:ga$window[2]
  Rtru <- mean(rt[[which]][rt$day %in% win], na.rm = TRUE)
  list(estimate_R = ga$R, truth_R = Rtru, error = ga$R - Rtru,
       in_ci = Rtru >= ga$R_ci[1] && Rtru <= ga$R_ci[2])
}

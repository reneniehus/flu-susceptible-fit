# phase1_cfr.R
#
# PHASE 1 -- early local exponential growth.
# QUESTION: "How deadly is it?"
#
# THE LEAN METHOD (Nishiura et al. 2009). The naive case fatality ratio, cumulative deaths / cumulative
# cases, is biased DOWN early in an epidemic because recent cases have not yet had time to die. The
# delay-adjusted (or "confirmed") CFR fixes this by shrinking the denominator to the cases whose
# outcome is already KNOWN in expectation, using the onset->death delay distribution f (CDF F):
#     known-outcome cases by day t  =  sum_{onset i <= t} cases_i * F(t - i)
#     cCFR(t)                       =  cumulative deaths(t) / known-outcome cases(t)
# As the epidemic matures F(t - i) -> 1 for most cases and the adjustment vanishes.
#
# WHAT IT DOES AND DOES NOT CORRECT. It corrects the DELAY bias only. It does NOT correct case
# UNDER-ASCERTAINMENT: if only a fraction rho of infections are detected as cases, the cCFR estimates
# deaths-per-DETECTED-case, which overstates the infection fatality ratio by ~ death_detection / rho.
# So the cCFR converges to the confirmed CFR, not the IFR -- the IFR needs seroprevalence (or the
# known ascertainment) on top. The playground makes both biases visible by scoring against truth.
#
# References
#   Nishiura H, Klinkenberg D, Roberts M, Heesterbeek JAP. Early epidemiological assessment of the
#     virulence of emerging infectious diseases in a population. PLoS ONE. 2009;4(8):e6852.
#   {cfr} R package (Epiverse-TRACE): cfr_static / cfr_rolling, which this emulates.
#
# Requires: analysis_common.R; R/epidist.R (discretise) for the onset->death delay.

# ---- |-expected fraction of cases-to-date whose outcome (death or survival) is already known ----
# cases : cases by onset day (ordered). o2d_cdf : onset->death CDF at lags 0,1,2,... as_of : cutoff day
# index. Returns u in [0,1]; the denominator of the adjusted CFR is u * cumulative cases.
known_outcome_fraction <- function(cases, day, o2d_cdf, as_of) {
  keep <- day <= as_of
  c_i  <- cases[keep]; lag <- as_of - day[keep]
  F_lag <- o2d_cdf[pmin(lag + 1, length(o2d_cdf))]        # F(as_of - onset), clamped to the CDF's support
  denom_known <- sum(c_i * F_lag)                          # cases with an expected-known outcome
  total <- sum(c_i)
  list(u = denom_known / max(total, 1), known = denom_known, total = total)
}

# ---- |-delay-adjusted (confirmed) CFR at a single cutoff, with a binomial CI ----
# input : as_analysis_input(sim). location : which location. as_of : the day to estimate at. series :
# which case series to use ("cases_by_onset" recommended -- it pairs with the onset->death delay).
cfr_static <- function(input, location, as_of = input$as_of,
                       series = c("cases_by_onset", "cases_by_report"), o2d = NULL, level = 0.95) {
  series  <- match.arg(series)
  ca      <- loc_series(input[[series]], location)
  de      <- loc_series(input$deaths, location)
  o2d_cdf <- cumsum(discretise(o2d %||% input$delays$onset_to_death))

  deaths_cum <- sum(de$deaths_by_date[de$day <= as_of])
  ko <- known_outcome_fraction(ca$cases, ca$day, o2d_cdf, as_of)

  naive <- deaths_cum / max(ko$total, 1)                  # cumulative deaths / cumulative cases
  cfr   <- deaths_cum / max(ko$known, 1)                  # delay-adjusted
  ci    <- .binom_ci(deaths_cum, ko$known, level)         # binomial CI on the adjusted ratio
  list(location = location, as_of = as_of, deaths = deaths_cum, cases = ko$total,
       known_outcome = ko$known, u = ko$u,
       cfr_naive = naive, cfr_adjusted = cfr, cfr_lower = ci[1], cfr_upper = ci[2])
}

# ---- |-delay-adjusted CFR as an expanding time series (naive vs adjusted, day by day) ----
# Shows the two curves converging: the naive CFR climbs from below as deaths catch up, the adjusted
# CFR is roughly flat from early on. cutoffs : the days to evaluate at (default every 7 days).
cfr_rolling <- function(input, location, cutoffs = NULL,
                        series = c("cases_by_onset", "cases_by_report"), o2d = NULL, level = 0.95) {
  series <- match.arg(series)
  ca <- loc_series(input[[series]], location)
  if (is.null(cutoffs)) cutoffs <- seq(min(ca$day) + 14, input$as_of, by = 7)
  rows <- lapply(cutoffs, function(t) {
    s <- cfr_static(input, location, as_of = t, series = series, o2d = o2d, level = level)
    data.frame(as_of = t, cfr_naive = s$cfr_naive, cfr_adjusted = s$cfr_adjusted,
               cfr_lower = s$cfr_lower, cfr_upper = s$cfr_upper, deaths = s$deaths, cases = s$cases)
  })
  do.call(rbind, rows)
}

# ---- |-Wilson score confidence interval for a ratio deaths / known ----
# The Wilson score interval, not the Wald p +/- z*se -- Wald under-covers and pins the lower bound at 0
# in exactly the small-count / small-p regime the early-epidemic CFR lives in (Brown, Cai & DasGupta
# 2001, Stat Sci 16:101-133).
.binom_ci <- function(deaths, known, level = 0.95) {
  if (known <= 0) return(c(NA_real_, NA_real_))
  p <- min(deaths / known, 1); n <- known; z <- stats::qnorm(1 - (1 - level) / 2)
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(max(centre - half, 0), min(centre + half, 1))
}

# ---- |-score the CFR against the true confirmed CFR and the true IFR ----
# The delay-adjusted CFR targets the CONFIRMED CFR: eventual detected deaths / eventual detected cases
# for the cohort. It should hit that; it should NOT hit the (much smaller) IFR, because case
# under-ascertainment is uncorrected -- and the gap to the IFR equals the under-ascertainment.
score_cfr <- function(sim, cfr, input = as_analysis_input(sim)) {
  loc <- cfr$location; as_of <- cfr$as_of
  # confirmed-CFR truth: eventual detected deaths-by-onset / eventual detected cases-by-onset, cohort <= as_of
  de <- loc_series(sim$observed$deaths, loc)
  ca <- loc_series(sim$truth$cases_by_onset, loc)                 # eventual detected cases by onset
  deaths_by_onset_cum <- sum(de$deaths_by_onset[de$day <= as_of])
  cases_eventual_cum  <- sum(ca$cases[ca$day <= as_of])
  confirmed_cfr <- deaths_by_onset_cum / max(cases_eventual_cum, 1)

  list(cfr_adjusted = cfr$cfr_adjusted, cfr_naive = cfr$cfr_naive,
       true_confirmed_cfr = confirmed_cfr, true_ifr = truth_ifr(sim),
       err_vs_confirmed = cfr$cfr_adjusted - confirmed_cfr,
       ratio_cfr_to_ifr = cfr$cfr_adjusted / truth_ifr(sim))     # the uncorrected ascertainment gap
}

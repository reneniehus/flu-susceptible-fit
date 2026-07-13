# phase1_rt.R
#
# PHASE 1 -- early local exponential growth (and beyond).
# QUESTION: "Is transmission growing or shrinking here, right now?"
#
# THE LEAN METHOD (Cori et al. 2013 -- the EpiEstim algorithm, self-contained). Model incidence as a
# renewal process: today's cases are Poisson with mean R_t * Lambda_t, where the total infectiousness
# Lambda_t = sum_{s>=1} w_s I_{t-s} convolves past incidence with the generation interval w. With a
# Gamma(a, b) prior on R, the posterior over a sliding window [t-tau+1, t] is conjugate:
#     R_t | data  ~  Gamma( a + sum_k I_k ,  rate = 1/b + sum_k Lambda_k )      (k in the window)
# giving a closed-form posterior mean and credible interval -- no MCMC, no dependencies. This is the
# leanest defensible Rt estimator and the one EpiEstim implements; EpiNow2 is the fuller (delay- and
# truncation-aware) alternative and a documented extension front.
#
# CAVEATS THE PLAYGROUND MAKES VISIBLE. (1) Feeding it CASES rather than infections assumes constant
# ascertainment -- a scaling-up testing regime inflates apparent R (the same ramp that biased Phase 0
# growth). (2) The most recent days are right-truncated (reports still arriving), biasing recent R
# down. (3) It should be fed cases by ONSET (or infection) date; by-report date smears the timing.
#
# References
#   Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and software to estimate time-varying
#     reproduction numbers during epidemics. Am J Epidemiol. 2013;178(9):1505-1512.
#
# Requires: analysis_common.R; R/epidist.R (discretise) for the generation interval.

# ---- |-Cori sliding-window R_t from an incidence series + generation interval ----
# incidence, day : aligned, CONTIGUOUS daily counts (cases by onset, or infections). gi_pmf :
# discretised generation interval (gi_pmf[1] = P(0) = 0). window : sliding-window length (days).
# prior_mean/prior_sd : Gamma prior on R (Cori default mean 5, sd 5). Returns a per-day data frame of
# the posterior mean and credible interval, attributed to the END of each window.
estimate_rt_cori <- function(incidence, day, gi_pmf, window = 7,
                             prior_mean = 5, prior_sd = 5, level = 0.95, min_signal = 5) {
  n <- length(incidence)
  # total infectiousness Lambda_t = sum_{s>=1} w_s I_{t-s}
  Lambda <- .total_infectiousness(incidence, gi_pmf)
  a0 <- (prior_mean / prior_sd)^2                        # Gamma prior shape
  b0 <- prior_sd^2 / prior_mean                          # Gamma prior scale (mean = a0*b0)

  lo <- (1 - level) / 2; hi <- 1 - lo
  out <- lapply(seq_len(n), function(t) {
    if (t <= window) return(NULL)                        # need a full window of past infectiousness
    idx  <- (t - window + 1):t
    sumI <- sum(incidence[idx]); sumL <- sum(Lambda[idx])
    # LOW-SIGNAL GUARD: with too few cases the Gamma posterior just returns the prior mean (~5 by
    # default), which looks like a confident estimate but is not. Refuse instead -- matching the
    # growth / variant tools, which stop() on too little signal. Critical for sub-critical / mild /
    # burnt-out epidemics where the sparse tail would otherwise read R ~ prior_mean on ~zero cases.
    if (sumL <= 0 || sumI < min_signal) return(NULL)
    shape <- a0 + sumI; rate <- 1 / b0 + sumL            # posterior Gamma(shape, rate)
    data.frame(day = day[t],
               Rt = shape / rate,
               Rt_lower = stats::qgamma(lo, shape, rate),
               Rt_upper = stats::qgamma(hi, shape, rate))
  })
  res <- do.call(rbind, out)
  # always return a typed data frame (0-row if nothing was estimable), never a bare NULL, so callers
  # and scoring fail cleanly rather than silently building a malformed object.
  if (is.null(res)) data.frame(day = integer(0), Rt = numeric(0),
                               Rt_lower = numeric(0), Rt_upper = numeric(0)) else res
}

# ---- |-total infectiousness Lambda_t = sum_{s>=1} w_s I_{t-s} ----
.total_infectiousness <- function(incidence, gi_pmf) {
  n <- length(incidence); D <- length(gi_pmf) - 1
  Lambda <- numeric(n)
  for (s in 1:D) if (gi_pmf[s + 1] > 0 && s < n)
    Lambda[(s + 1):n] <- Lambda[(s + 1):n] + gi_pmf[s + 1] * incidence[1:(n - s)]
  Lambda
}

# ---- |-Phase-1 Rt analysis for one location from the observed onset-date cases ----
# input : as_analysis_input(sim). Uses cases by onset by default (the right series for a GI-based
# renewal estimate). `truncate` drops the last few days, where right-truncation biases R down.
# `window = NULL` DERIVES the sliding window from the generation interval (~1.5 x GI mean), so it is
# ~8 days for COVID but ~4 for fast flu and ~18 for measles -- a fixed 7-day window is 2.7 generations
# for flu (over-smoothed, lags the turnover) and 0.6 for measles (unstable). Pass a number to override.
rt_analysis <- function(input, location, window = NULL, gi = NULL,
                        series = c("cases_by_onset", "cases_by_report"),
                        prior_mean = 5, prior_sd = 5, truncate = 0, min_signal = 5) {
  series <- match.arg(series)
  gi_dist <- gi %||% input$delays$generation_interval
  gi_pmf  <- discretise(gi_dist, boundary = "cori")
  if (is.null(window)) window <- max(3L, round(1.5 * epidist_mean(gi_dist)))   # GI-scaled default
  d  <- loc_series(input[[series]], location)
  if (truncate > 0) d <- d[d$day <= (max(d$day) - truncate), ]
  est <- estimate_rt_cori(d$cases, d$day, gi_pmf, window = window,
                          prior_mean = prior_mean, prior_sd = prior_sd, min_signal = min_signal)
  est$location <- location
  est
}

# ---- |-score estimated R_t against the true REALIZED instantaneous R_t (aligned by day) ----
# The right target is `Rt_effective` -- the strain-mix- and depletion-weighted reproduction number the
# epidemic actually ran at -- not the nominal wild-type schedule (which ignores the variant and
# depletion). Step changes are smeared by the sliding window, so agreement is on level and direction,
# not sharp edges. `which = "Rt"` scores against the nominal schedule instead, if wanted.
score_rt <- function(sim, rt_est, which = c("Rt_effective", "Rt")) {
  which <- match.arg(which)
  if (!nrow(rt_est)) return(list(n = 0, note = "no Rt was estimable (too little signal)"))
  loc   <- rt_est$location[1]
  truth <- truth_rt(sim, loc)
  Rtrue <- truth[[which]][match(rt_est$day, truth$day)]
  s <- pp_score(rt_est$Rt, Rtrue, rt_est$Rt_lower, rt_est$Rt_upper)
  s$location <- loc
  s
}

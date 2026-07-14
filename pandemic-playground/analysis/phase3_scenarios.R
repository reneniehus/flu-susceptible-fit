# phase3_scenarios.R
#
# PHASE 3 -- sustained transmission, later waves, endemic-ish.
# QUESTION: "Given waning immunity / a new variant / a booster campaign, what might next season look
#            like, and does boosting help?"
#
# THE LEAN METHOD. A parsimonious SIRS model, integrated with a plain RK4 stepper (no deSolve/odin
# needed), run over a small FACTORIAL of the three levers that actually move a next season:
#   - waning rate omega   (how fast immunity is lost, 1/duration)
#   - transmissibility R0 (a new, fitter variant raises it)
#   - booster uptake      (a campaign that moves a fraction of susceptibles into the immune class)
# We report each scenario's peak prevalence, peak timing and cumulative incidence RELATIVE to a reference
# scenario -- differences, not fragile absolute levels -- and ENSEMBLE across the most uncertain axis
# (waning) so a headline like "boosting cuts the peak by ~X%" carries its assumption spread with it.
# This mirrors a scenario-modelling-hub target: a tidy table of scenarios and comparable outcomes.
#
# WHY RELATIVE, ENSEMBLED. Absolute next-season numbers depend on parameters no one knows in advance
# (see PROJECT_SCOPE of the parent flu project); relative differences between scenarios that share
# those unknowns are far more robust, and averaging over the waning assumption states the uncertainty
# instead of hiding it. Seed the initial immune fraction from where the current epidemic actually
# left the population (from a simulated -- or real -- run) to ground "next season" in reality.
#
# THIS IS THE LEAST WELL-POSED TOOL IN THE BOX -- BY CONSTRUCTION, NOT BY SLOPPINESS. A long-horizon SIRS
# projection ("does boosting help next season?") sits furthest outside the envelope of anything we have
# observed: over a season, behaviour, a new variant's properties, and waning are all things the model
# ASSUMES rather than measures, and small differences in those assumptions swamp the mechanistic detail.
# That is precisely why the output is engineered to be modest -- RELATIVE outcomes, ENSEMBLED over the
# most uncertain axis -- and why an absolute "peak of X admissions in week W next winter" from this tool
# should be treated as illustrative arithmetic, not a prediction. The compartmental model is confidently
# wrong OUTSIDE its envelope for the same reason a data-driven / foundation model would be confidently
# wrong here too (it has no analogous season to lean on). Use this to RANK levers ("boosting the elderly
# beats a broad campaign under most waning assumptions"), never to quote a number. (See reflections.md:
# "which idealisations are well-posed", and why honest out-of-envelope detection is the real prize.)
#
# References
#   Keeling MJ, Rohani P. Modeling Infectious Diseases in Humans and Animals. Princeton; 2008. (SIRS)
#   ECDC RespiCompass / US Scenario Modeling Hub -- the scenario-grid, relative-outcome format.
#
# Requires: analysis_common.R.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### A lean SIRS integrator (base-R RK4) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-integrate an SIRS model in proportions with RK4 ----
# dS/dt = -beta S I + omega R ; dI/dt = beta S I - gamma I ; dR/dt = gamma I - omega R ; and a
# cumulative-incidence accumulator dC/dt = beta S I (infection events per capita; equals the final
# size 1 - S_end for a pure SIR, but can exceed 1 under waning). beta = R0 * gamma; S+I+R = 1.
# Integrates with a fixed number of RK4 sub-steps PER DAY, so every integer day lands exactly and the
# returned trajectory always spans day 0..days (no dropped final day for a non-divisor dt).
sirs_integrate <- function(R0, gamma, omega, S0, I0, days, dt = 0.25) {
  beta <- R0 * gamma
  deriv <- function(y) {
    S <- y[1]; I <- y[2]; R <- y[3]; inc <- beta * S * I
    c(-inc + omega * R, inc - gamma * I, gamma * I - omega * R, inc)
  }
  steps_per_day <- max(1L, as.integer(round(1 / dt)))    # snap dt to an integer number of sub-steps/day
  h <- 1 / steps_per_day
  out <- matrix(NA_real_, nrow = days + 1, ncol = 4)
  y <- c(S0, I0, 1 - S0 - I0, 0); out[1, ] <- y
  for (d in seq_len(days)) {
    for (s in seq_len(steps_per_day)) {
      k1 <- deriv(y); k2 <- deriv(y + h / 2 * k1)
      k3 <- deriv(y + h / 2 * k2); k4 <- deriv(y + h * k3)
      y  <- pmax(y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 0)   # step + numerical guard
    }
    out[d + 1, ] <- y
  }
  data.frame(day = 0:days, S = out[, 1], I = out[, 2], R = out[, 3], C = out[, 4])
}

# ---- |-summary outcomes of one SIRS trajectory ----
# cumulative_incidence is infection EVENTS per capita over the horizon (the C accumulator). It equals
# the final size 1 - S_end for a pure SIR, but CAN EXCEED 1 under waning immunity, because a person
# can be reinfected within the season -- so it is a burden measure, not a "fraction ever infected".
.sirs_outcomes <- function(traj) {
  list(peak_prevalence = max(traj$I), peak_day = traj$day[which.max(traj$I)],
       cumulative_incidence = utils::tail(traj$C, 1))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### The scenario factorial ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-run a factorial of SIRS scenarios and report outcomes relative to a reference ----
# grid : a list of the three levers, each a named numeric vector, e.g.
#   list(waning = c(slow = 1/365, fast = 1/120),
#        R0 = c(wildtype = 1.3, variant = 1.8),
#        uptake = c(none = 0, campaign = 0.4))
# immune0 : starting immune fraction (seed from a run's end-of-wave susceptibles: 1 - S_end).
# gamma sets the recovery rate and thus the WHOLE epidemic speed; pass 1/mean(generation interval) for
# the pathogen (default 1/5 is COVID-like -- a fast flu wants ~1/3, measles ~1/12). I0, days : constants.
# reference : the scenario names to treat as the comparison baseline.
sirs_scenarios <- function(grid, immune0, gamma = 1 / 5, I0 = 1e-4, days = 365,
                           reference = c(waning = names(grid$waning)[1],
                                         R0 = names(grid$R0)[1], uptake = names(grid$uptake)[1])) {
  combos <- expand.grid(waning = names(grid$waning), R0 = names(grid$R0), uptake = names(grid$uptake),
                        stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    w <- combos$waning[i]; v <- combos$R0[i]; u <- combos$uptake[i]
    uptake <- grid$uptake[[u]]
    S0 <- (1 - immune0) * (1 - uptake)                  # a booster campaign protects `uptake` of the susceptibles
    traj <- sirs_integrate(R0 = grid$R0[[v]], gamma = gamma, omega = grid$waning[[w]],
                           S0 = S0, I0 = I0, days = days)
    o <- .sirs_outcomes(traj)
    data.frame(scenario = sprintf("%s | %s | %s", w, v, u), waning = w, R0 = v, uptake = u,
               peak_prevalence = o$peak_prevalence, peak_day = o$peak_day,
               cumulative_incidence = o$cumulative_incidence)
  })
  res <- do.call(rbind, rows)

  ref_id <- sprintf("%s | %s | %s", reference["waning"], reference["R0"], reference["uptake"])
  ref <- res[res$scenario == ref_id, ]
  res$rel_peak       <- res$peak_prevalence / ref$peak_prevalence          # vs the reference scenario
  res$rel_incidence  <- res$cumulative_incidence / ref$cumulative_incidence
  res
}

# ---- |-ensemble over the waning axis: average each (R0, uptake) cell across waning assumptions ----
# Turns the factorial into an assumption-averaged headline with a spread, e.g. "the variant + no
# booster raises the peak ~X% (range across waning) vs wildtype + booster".
sirs_ensemble <- function(res) {
  agg <- stats::aggregate(cbind(rel_peak, rel_incidence) ~ R0 + uptake, data = res,
                          FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
  # aggregate() with a vector-returning FUN makes agg$rel_peak a MATRIX (columns mean/min/max), so
  # agg$rel_peak[, "mean"] below is a sub-column selection, not a typo
  data.frame(R0 = agg$R0, uptake = agg$uptake,
             rel_peak_mean = agg$rel_peak[, "mean"], rel_peak_min = agg$rel_peak[, "min"],
             rel_peak_max = agg$rel_peak[, "max"],
             rel_incidence_mean = agg$rel_incidence[, "mean"])
}

# ---- |-"does boosting help?" -- the campaign-vs-no-campaign effect for a given variant ----
# Returns the relative peak / cumulative-incidence reduction from the booster campaign, ensembled
# across waning.
boosting_effect <- function(res, variant_level, none = "none", campaign = NULL) {
  if (is.null(campaign)) campaign <- setdiff(unique(res$uptake), none)[1]
  a <- res[res$R0 == variant_level & res$uptake == campaign, ]
  b <- res[res$R0 == variant_level & res$uptake == none, ]
  m <- merge(a, b, by = "waning", suffixes = c("_boost", "_none"))
  list(variant = variant_level, campaign = campaign,
       peak_reduction      = 1 - mean(m$peak_prevalence_boost / m$peak_prevalence_none),
       incidence_reduction = 1 - mean(m$cumulative_incidence_boost / m$cumulative_incidence_none))
}

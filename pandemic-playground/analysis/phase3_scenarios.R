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
# We report each scenario's peak prevalence, peak timing and attack rate RELATIVE to a reference
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
# References
#   Keeling MJ, Rohani P. Modeling Infectious Diseases in Humans and Animals. Princeton; 2008. (SIRS)
#   ECDC RespiCompass / US Scenario Modeling Hub -- the scenario-grid, relative-outcome format.
#
# Requires: analysis_common.R.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### A lean SIRS integrator (base-R RK4) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-integrate an SIRS model in proportions with RK4 ----
# dS/dt = -beta S I + omega R ; dI/dt = beta S I - gamma I ; dR/dt = gamma I - omega R.
# beta = R0 * gamma. State is (S, I, R), summing to 1. Returns day-indexed S, I, R trajectories.
sirs_integrate <- function(R0, gamma, omega, S0, I0, days, dt = 0.25) {
  beta <- R0 * gamma
  deriv <- function(y) {
    S <- y[1]; I <- y[2]; R <- y[3]
    c(-beta * S * I + omega * R,
       beta * S * I - gamma * I,
       gamma * I - omega * R)
  }
  n_steps <- as.integer(days / dt)
  y <- c(S0, I0, 1 - S0 - I0)
  out <- matrix(NA_real_, nrow = days + 1, ncol = 3); out[1, ] <- y; keep <- 1
  next_keep_t <- 1
  for (k in seq_len(n_steps)) {
    k1 <- deriv(y); k2 <- deriv(y + dt / 2 * k1)
    k3 <- deriv(y + dt / 2 * k2); k4 <- deriv(y + dt * k3)
    y  <- y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    y  <- pmax(y, 0)                                    # numerical guard
    if ((k * dt) >= next_keep_t) { keep <- keep + 1; out[keep, ] <- y; next_keep_t <- next_keep_t + 1 }
  }
  data.frame(day = 0:(keep - 1), S = out[1:keep, 1], I = out[1:keep, 2], R = out[1:keep, 3])
}

# ---- |-summary outcomes of one SIRS trajectory ----
.sirs_outcomes <- function(traj) {
  list(peak_prevalence = max(traj$I), peak_day = traj$day[which.max(traj$I)],
       attack_rate = sum(traj$I))                       # sum of I over days ~ cumulative infection-days (burden proxy)
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
# gamma, I0, days : SIRS constants. reference : the scenario names to treat as the comparison baseline.
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
               peak_prevalence = o$peak_prevalence, peak_day = o$peak_day, attack_rate = o$attack_rate)
  })
  res <- do.call(rbind, rows)

  ref_id <- sprintf("%s | %s | %s", reference["waning"], reference["R0"], reference["uptake"])
  ref <- res[res$scenario == ref_id, ]
  res$rel_peak   <- res$peak_prevalence / ref$peak_prevalence      # vs the reference scenario
  res$rel_attack <- res$attack_rate / ref$attack_rate
  res
}

# ---- |-ensemble over the waning axis: average each (R0, uptake) cell across waning assumptions ----
# Turns the factorial into an assumption-averaged headline with a spread, e.g. "the variant + no
# booster raises the peak ~X% (range across waning) vs wildtype + booster".
sirs_ensemble <- function(res) {
  agg <- stats::aggregate(cbind(rel_peak, rel_attack) ~ R0 + uptake, data = res,
                          FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
  data.frame(R0 = agg$R0, uptake = agg$uptake,
             rel_peak_mean = agg$rel_peak[, "mean"], rel_peak_min = agg$rel_peak[, "min"],
             rel_peak_max = agg$rel_peak[, "max"],
             rel_attack_mean = agg$rel_attack[, "mean"])
}

# ---- |-"does boosting help?" -- the campaign-vs-no-campaign effect for a given variant ----
# Returns the relative peak/attack reduction from the booster campaign, ensembled across waning.
boosting_effect <- function(res, variant_level, none = "none", campaign = NULL) {
  if (is.null(campaign)) campaign <- setdiff(unique(res$uptake), none)[1]
  a <- res[res$R0 == variant_level & res$uptake == campaign, ]
  b <- res[res$R0 == variant_level & res$uptake == none, ]
  m <- merge(a, b, by = "waning", suffixes = c("_boost", "_none"))
  list(variant = variant_level, campaign = campaign,
       peak_reduction   = 1 - mean(m$peak_prevalence_boost / m$peak_prevalence_none),
       attack_reduction = 1 - mean(m$attack_rate_boost / m$attack_rate_none))
}

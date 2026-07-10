# phase0_catchment.R
#
# PHASE 0 -- before local introduction, at the source X.
# QUESTION: "How much are we missing at source -- how big is it really?"
#
# THE LEAN METHOD. Exported cases detected abroad are a window onto the source epidemic that does NOT
# depend on the source's own (poor, early) surveillance. If a fraction of the source population flies
# to a set of well-surveilled catchment countries each day, and those countries detect imported
# infections with known probability, then the number of infectious people at the source is
#     prevalence_X  ~  (detected exports to the catchment) / (sum_c travel-fraction_c * detection_c)
# where travel-fraction_c = volume_{X->c} / N_X and detection_c is c's surveillance quality. Anchoring
# on WELL-SURVEILLED destinations (high detection, which we treat as known) makes the denominator
# trustworthy. Comparing the back-calculated size to the source's OWN reported cases reveals the
# under-ascertainment multiple -- "how much bigger than it looks".
#
# This is the arithmetic version; a small Poisson model (detected ~ Poisson(volume * detection * prev))
# is the natural next step and is a documented extension front.
#
# References
#   Fraser C, et al. Pandemic potential of a strain of influenza A (H1N1). Science. 2009;324:1557-1561.
#   De Salazar PM, et al. Using predicted imports of 2019-nCoV cases to determine locations that may
#     not be identifying all imported cases. medRxiv/Emerg Infect Dis. 2020.  (catchment logic)
#
# Requires: analysis_common.R.

# ---- |-back-calculate source infectious prevalence from detected exports + travel + detection ----
# input : as_analysis_input(sim). window : day indices over which to aggregate (an early window, while
# the source epidemic still dominates). min_surveillance : anchor only on countries at least this well
# surveilled (their detection is trusted). source_pop : X's population (needed to turn a fraction into
# a count; defaults to the sim's known value if a sim is passed via `source_pop`).
catchment_backcalc <- function(input, window, min_surveillance = 0.8, source_pop) {
  imp <- input$detected_imports
  vol <- input$flight_volumes
  # anchor countries: those whose surveillance quality clears the threshold (trusted detection)
  sq  <- tapply(imp$surveillance_quality, imp$country, function(x) x[1])
  anchor <- names(sq)[sq >= min_surveillance]
  if (!length(anchor)) stop("catchment_backcalc: no countries meet min_surveillance; lower the threshold")

  in_win <- function(d) d$day %in% window & d$country %in% anchor
  imp_w  <- imp[in_win(imp), ]
  vol_w  <- vol[in_win(vol), ]

  detected_exports <- sum(imp_w$detected_imports)
  # sum over anchor countries & window of (volume * detection): the effective "sampling effort"
  sampling <- vol_w |>
    dplyr::inner_join(dplyr::distinct(imp[, c("country", "surveillance_quality")]), by = "country") |>
    dplyr::mutate(effort = volume * surveillance_quality) |>
    dplyr::summarise(x = sum(effort)) |> dplyr::pull(x)

  est_prev_frac  <- detected_exports / sampling           # volume-weighted average source prevalence fraction
  est_prevalence <- est_prev_frac * source_pop            # infectious people at X (count)

  list(window = range(window), n_anchor = length(anchor), anchor = anchor,
       detected_exports = detected_exports,
       est_prevalence_frac = est_prev_frac, est_prevalence = est_prevalence,
       source_pop = source_pop)
}

# ---- |-score the catchment estimate against the true source prevalence over the window ----
# Truth: the true source infectious prevalence, volume-weighted-averaged over the same window (the
# quantity the catchment sum actually estimates -- a STOCK of infectious people). Also reports the
# TRUE source under-ascertainment (a clean flow-vs-flow ratio: true new infections vs reported cases
# over the window) -- the reality the exported-case method exists to reveal, read straight from truth.
score_catchment <- function(sim, cb, input = as_analysis_input(sim)) {
  win  <- cb$window[1]:cb$window[2]
  prev <- sim$truth$prevalence
  vol  <- sim$latent$flights$true                          # true volumes (for the exact weighting)
  w    <- rowSums(vol)[win + 1]                             # total outbound volume per day (weights)
  true_prev_frac <- stats::weighted.mean(prev$prevalence_frac[win + 1], w)
  true_prev      <- true_prev_frac * cb$source_pop

  # the reality: true source infections vs the source's own reported cases over the window (flow/flow)
  true_inf  <- truth_infections(sim, input$source_code)
  reported  <- loc_series(input$cases_by_report, input$source_code)
  true_infections_win <- sum(true_inf$infections[true_inf$day %in% win])
  reported_win        <- sum(reported$cases[reported$day %in% win])

  list(est_prevalence = cb$est_prevalence, true_prevalence = true_prev,
       ratio = cb$est_prevalence / true_prev,
       reported_source_cases = reported_win, true_source_infections = true_infections_win,
       true_underascertainment = true_infections_win / max(reported_win, 1))
}

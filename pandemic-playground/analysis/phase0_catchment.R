# phase0_catchment.R
#
# PHASE 0 -- before local introduction, at the source X.
# QUESTION: "How much are we missing at source -- how big is it really?"
#
# THE METHOD (travellers as a sentinel sample of the source). Exported cases detected abroad are a
# window onto the source epidemic that does NOT depend on the source's own (poor, early) surveillance.
# It is a simple balance equation: if a fraction of the source population flies to a set of
# well-surveilled catchment countries each day, and those countries detect imported infections, then
#     prevalence_X  ~  (detected exports to the catchment) / (sum_c travel-fraction_c * detection_c)
# with travel-fraction_c = volume_{X->c} / N_X. Anchoring on WELL-SURVEILLED destinations makes the
# denominator trustworthy; comparing the back-calculated size to the source's OWN reported cases gives
# the under-ascertainment multiple -- "how much bigger than it looks".
#
# THIS WORKED, TWICE, ON REAL EMERGENCIES. In Wuhan (Jan 2020) three cases exported to Thailand/Japan
# plus outbound air volume back-calculated ~1,700 cases when ~40 were officially reported -- a ~40x
# under-count (Imperial College; Wu, Leung & Leung then fitted ~75,000 infections). For 2009 H1N1 the
# same logic, run from Mexico's exported cases, sized the outbreak and anchored the CFR and R0 (Fraser
# et al.). The honest output was always a RANGE, because two corrections must be made explicit:
#   (1) REPRESENTATIVENESS (De Salazar et al. 2020). Short-stay VISITORS under-sample source prevalence
#       relative to residents by a factor V = 1 - exp(-(r + gamma) * d), with length of stay d, source
#       growth rate r and recovery rate gamma. Long-stay travellers give V ~ 1; a 3-day visit ~ 0.5.
#       Using short-stay travellers biases the source estimate DOWNWARD -- divide by V to correct.
#   (2) DIFFERENTIAL DETECTION (Niehus et al. 2020). Benchmark every destination against the best
#       surveillance; global detection was only ~38% of the best, so raw exported-case counts are too
#       low. A `detection_scale` < 1 (the anchors' detection relative to a perfect system) scales up.
# Get these wrong and the source size is biased -- almost always downward -- so report the range with
# the assumptions stated. (In THIS simulator travellers are representative and detection = surveillance
# quality exactly, so the uncorrected estimate is unbiased; the corrections are the real-data levers,
# shown here as an explicit sensitivity range via catchment_range().)
#
# References
#   Imperial College COVID-19 Response Team. Report 2: estimating the potential total number of novel
#     coronavirus cases in Wuhan City from those exported. 2020.
#   Wu JT, Leung K, Leung GM. Nowcasting and forecasting the potential domestic and international spread
#     of 2019-nCoV. Lancet. 2020;395:689-697.
#   Fraser C, et al. Pandemic potential of a strain of influenza A (H1N1). Science. 2009;324:1557-1561.
#   De Salazar PM, et al. Using predicted imports of 2019-nCoV cases... Emerg Infect Dis. 2020.
#     (representativeness V)   Niehus R, et al. Lancet Digit Health / eLife 2020. (differential detection)
#
# Requires: analysis_common.R.

# ---- |-visitor-vs-resident representativeness factor V = 1 - exp(-(r + gamma) * d) ----
# The ratio of prevalence in travellers (mean stay d) to residents. Long stays -> V ~ 1 (representative);
# short visits under-sample the growing prevalence. Divide a source estimate by V to correct upward.
visitor_resident_factor <- function(stay_days, growth_rate, recovery_rate) {
  1 - exp(-(growth_rate + recovery_rate) * stay_days)
}

# ---- |-back-calculate source infectious prevalence from detected exports + travel + detection ----
# input : as_analysis_input(sim). window : day indices to aggregate over (an early window, while the
# source epidemic dominates). min_surveillance : anchor only on countries at least this well surveilled.
# source_pop : X's population. stay_days/growth_rate/recovery_rate : optional -- apply the De Salazar
# representativeness correction V (NULL -> travellers assumed representative, V = 1). detection_scale :
# the anchors' detection relative to a perfect system (Niehus; 1 -> take surveillance quality as absolute).
catchment_backcalc <- function(input, window, min_surveillance = 0.8, source_pop,
                               stay_days = NULL, growth_rate = NULL, recovery_rate = NULL,
                               detection_scale = 1) {
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

  # apply the two real-data corrections: representativeness V and the detection benchmark. Both sit in
  # the denominator, so a smaller V or detection_scale raises the estimate (the usual direction).
  V <- if (!is.null(stay_days)) visitor_resident_factor(stay_days, growth_rate, recovery_rate) else 1
  est_prev_frac  <- detected_exports / (sampling * detection_scale * V)   # source prevalence fraction
  est_prevalence <- est_prev_frac * source_pop                            # infectious people at X (count)

  list(window = range(window), n_anchor = length(anchor), anchor = anchor,
       detected_exports = detected_exports, V = V, detection_scale = detection_scale,
       est_prevalence_frac = est_prev_frac, est_prevalence = est_prevalence,
       source_pop = source_pop)
}

# ---- |-the source size as an explicit RANGE across the correction assumptions ----
# Sweeps a plausible span of traveller stay length and anchor detection to give low / central / high
# source-size estimates -- the honest output form (as the Wuhan estimates were presented). Both
# corrections only ever push the estimate UP (shorter stays and poorer detection mean each observed
# export stands for MORE at source), so the ordering is fixed:
#   low     = uncorrected  -- travellers representative (V=1), anchors detect everything (scale=1). Floor.
#   central = moderate stay + partial detection            -- the working point.
#   high    = shortest stay + poorest detection            -- ceiling.
# `growth_rate` and `recovery_rate` set the representativeness curve; supply the Phase-0 growth estimate.
catchment_range <- function(input, window, source_pop, growth_rate, recovery_rate,
                            min_surveillance = 0.8,
                            stay_days = c(long = 30, moderate = 7, short = 3),
                            detection_scale = c(high = 1, low = 0.5)) {
  size_at <- function(stay, det) catchment_backcalc(
    input, window, min_surveillance, source_pop, stay_days = stay,
    growth_rate = growth_rate, recovery_rate = recovery_rate, detection_scale = det)$est_prevalence

  low     <- catchment_backcalc(input, window, min_surveillance, source_pop)$est_prevalence  # V=1, det=1
  high    <- size_at(min(stay_days), min(detection_scale))                     # short stay + poor detection
  central <- size_at(stats::median(stay_days), mean(detection_scale))          # moderate stay + partial det
  list(central = central, low = low, high = high,
       assumptions = list(stay_days = stay_days, detection_scale = detection_scale,
                          growth_rate = growth_rate, recovery_rate = recovery_rate))
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
  # day indices are 0-based; +1 maps them to the 1-based rows of the volume / prevalence matrices
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

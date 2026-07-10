# simulate_flights.R
#
# Stage 2: the daily air-travel volume from the source X to each destination country. Volumes are
# built from a few interpretable pieces -- a baseline that scales with destination population, a
# persistent per-country connectivity multiplier (noise), a mild annual seasonality, and a travel-ban
# shock after a configurable date. This is the conduit that turns the source epidemic into
# importations, and (via its OBSERVED, noisy copy) the covariate the importation-risk analysis
# regresses on.
#
# We emit BOTH the TRUE volumes (what importation actually depends on) and an OBSERVED version blurred
# by multiplicative measurement noise -- so an analyst never sees the exact volumes, exactly as in
# reality. Keeping the two separate is the whole point of the truth/observation split.
#
# volume_{X->c, t} = baseline_scale * (pop_c / 1e6)^power * connectivity_c
#                    * seasonality(t) * ban(t)
#
# Requires: R/utils.R sourced first.

# ---- |-simulate true + observed daily flight volumes X -> each country ----
# cfg : config; par : draw_parameters(cfg). Uses the current RNG state (for the observation noise).
simulate_flights <- function(cfg, par) {
  n    <- cfg$n_days
  cc   <- cfg$countries
  fl   <- cfg$flights
  day  <- par$grid$day

  pop_m       <- par$country_pop / 1e6
  connectivity <- fl$country_noise[cc]                              # persistent per-country multiplier
  base_c      <- fl$baseline_scale * (pop_m ^ fl$population_power) * connectivity   # per-country baseline

  season <- seasonality_multiplier(day, fl$seasonality_amp, fl$seasonality_peak_day)  # length n

  ban <- rep(1, n)                                                 # travel-ban shock (multiplicative)
  if (!is.na(par$travel_ban_day)) ban[(par$travel_ban_day + 1):n] <- fl$travel_ban_factor

  # outer product day-modulation x per-country baseline -> n x n_countries TRUE volume matrix
  modulation <- season * ban                                       # length n (same for every country)
  true_volume <- outer(modulation, base_c)                         # n x n_countries
  colnames(true_volume) <- cc

  observed_volume <- matrix(rlnorm_measure(as.vector(true_volume), fl$observation_cv),
                            nrow = n, dimnames = list(NULL, cc))    # blurred copy the analyst sees

  list(true = true_volume, observed = observed_volume,
       baseline = base_c, connectivity = connectivity)
}

# draw_parameters.R
#
# Stage 0 of the pipeline: resolve a config into the concrete, ready-to-use quantities the rest of
# the simulator consumes. This is the ONE place the `epidist` delay distributions become discretised
# day-indexed PMFs, the age-structured IFR (if any) collapses to an effective IFR, and the calendar
# dates (variant introduction, travel ban) become integer day indices. Everything here is a
# DETERMINISTIC function of the config -- no random draws -- so it can be inspected and unit-tested on
# its own; all the stochasticity lives in the later stages under the run seed.
#
# Requires: R/epidist.R, R/utils.R sourced first.

# ---- |-resolve a config into discretised PMFs, an effective IFR, and day indices ----
# Returns a `par` list carried through the whole pipeline alongside the config.
draw_parameters <- function(cfg) {
  grid <- sim_days(cfg)

  # ---- delay PMFs (the generation interval uses the "no same-day" boundary; the rest keep day 0) ----
  pmf <- list(
    gi                 = discretise(cfg$delays$generation_interval, boundary = "cori"),
    incubation         = discretise(cfg$delays$incubation),
    onset_to_death     = discretise(cfg$delays$onset_to_death),
    onset_to_report    = discretise(cfg$delays$onset_to_report),
    onset_to_admission = if (!is.null(cfg$delays$onset_to_admission))
                           discretise(cfg$delays$onset_to_admission) else NULL
  )

  # ---- effective IFR: use the age-structured table if supplied, else the scalar ----
  # The DGP is not age-stratified (an extension front), so an age-structured IFR is reduced to a
  # single effective rate by the population age weights -- exact for a non-age-stratified epidemic.
  ifr_eff <- if (!is.null(cfg$ifr_age)) {
    a <- cfg$ifr_age
    stopifnot(all(c("ifr", "weight") %in% names(a)))
    sum(a$ifr * a$weight) / sum(a$weight)
  } else cfg$ifr

  # ---- number of strains and the calendar -> day-index conversions ----
  n_strains       <- if (isTRUE(cfg$variant$enabled)) 2L else 1L
  variant_day     <- if (n_strains == 2L) date_to_day(cfg$variant$intro_date, cfg) else NA_integer_
  travel_ban_day  <- if (!is.null(cfg$flights$travel_ban_date)) date_to_day(cfg$flights$travel_ban_date, cfg) else NA_integer_

  # ---- destination populations, aligned to the country vector ----
  geo  <- cfg$geography
  popn <- stats::setNames(geo$population_m[match(cfg$countries, geo$country)] * 1e6, cfg$countries)
  if (anyNA(popn)) stop("draw_parameters: some countries have no population in cfg$geography")

  list(
    grid            = grid,
    pmf             = pmf,
    ifr_eff         = ifr_eff,
    n_strains       = n_strains,
    variant_day     = variant_day,
    travel_ban_day  = travel_ban_day,
    country_pop     = popn,
    n_countries     = length(cfg$countries)
  )
}

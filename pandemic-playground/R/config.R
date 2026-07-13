# config.R
#
# The single knob-set for a playground run. `default_config()` returns a nested list holding EVERY
# assumption -- dates, the source location, the destination countries, the delay distributions, the
# reproduction-number trajectories, ascertainment, IFR, the flight model, per-country surveillance
# quality, and the variant -- so a whole synthetic pandemic is reproducible from one object plus a
# seed. Nothing downstream reads a hard-coded constant; it all comes from here. `validate_config()`
# checks the object before a run so mistakes fail loudly and early.
#
# DESIGN. The config is a plain list (like the parent repo's settings()), not an S4/R6 object, so it
# is trivial to read, print, tweak in a script, or serialise. Every field is documented inline with
# WHY it exists and what it drives. To change an experiment you copy `default_config()`, edit a
# field, and re-run -- see demo/ for worked examples. The delay distributions are `epidist` objects
# (R/epidist.R), so a curated {epiparameter} distribution drops straight in via as_epidist().
#
# Requires: R/epidist.R sourced first.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Reference geography: the EU/EEA destination set ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-EU/EEA countries with populations (millions, ~2023 Eurostat rounding) ----
# The default destination set. Population drives baseline flight connectivity and each country's
# susceptible pool, so it is the one piece of "real" geography the simulator needs. Kept as a small
# committed table -- no external data, fully reproducible. (27 EU + IS, NO, LI = the EEA.)
eu_eea_countries <- function() {
  data.frame(
    country = c("AT","BE","BG","HR","CY","CZ","DK","EE","FI","FR","DE","GR","HU","IE","IT",
                "LV","LT","LU","MT","NL","PL","PT","RO","SK","SI","ES","SE","IS","NO","LI"),
    population_m = c(9.10,11.75,6.45,3.85,0.92,10.83,5.94,1.37,5.56,68.17,84.36,10.41,9.60,5.19,58.85,
                     1.88,2.86,0.66,0.54,17.81,36.75,10.47,19.05,5.43,2.12,48.06,10.52,0.38,5.52,0.04),
    stringsAsFactors = FALSE)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Trajectory helpers: change-point step functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-a change-point (step) schedule for a time-varying quantity ----
# `t` are day indices (0-based, ascending, first MUST be 0); `value` is the level that holds from
# each change-point until the next. `step_schedule(t = c(0,40,70), value = c(2.4,1.2,0.9))` reads as
# "R = 2.4 from day 0, drops to 1.2 at day 40, to 0.9 at day 70". Evaluated by step_at() (R/utils.R).
# Reproduction numbers, ascertainment and travel shocks all use this one representation.
step_schedule <- function(t, value) {
  stopifnot(length(t) == length(value), t[1] == 0, !is.unsorted(t))
  list(t = as.numeric(t), value = as.numeric(value))
}

# ---- |-scale a reproduction-number schedule for the variant (x (1 + fitness)) ----
# The variant strain inherits a location's R trajectory times its transmissibility (fitness) advantage.
# Kept in one place so the "+fitness" rule is defined once and reused wherever a variant is set up (the
# source, the countries, and the realized-Rt truth).
scale_rt_for_variant <- function(base_rt, fitness) {
  base_rt$value <- base_rt$value * (1 + fitness)
  base_rt
}

# ---- |-build a per-country R_t schedule from a shared template + a country intervention ----
# Every country starts at `R_start`, then a control measure on `intervention_day` steps it to
# `R_post`. Countries differ only in WHEN they intervene (staggered response) -- the leanest way to
# get heterogeneous-but-comparable local epidemics. Override any single country by editing the list.
make_country_rt <- function(countries, R_start = 1.8, R_post = 0.85,
                            intervention_days) {
  stats::setNames(lapply(seq_along(countries), function(i)
    step_schedule(t = c(0, intervention_days[i]), value = c(R_start, R_post))), countries)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### The default configuration ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-a complete, reproducible default pandemic configuration ----
# Returns the nested list consumed by simulate_pandemic(). Every field is a lever; sensible
# respiratory-pathogen defaults are chosen so an out-of-the-box run produces a recognisable outbreak:
# a source epidemic at X that seeds staggered EU/EEA waves, with deaths, reporting delays, a travel
# ban and an optional fitter variant. Edit a copy of this to design an experiment.
default_config <- function() {
  geo   <- eu_eea_countries()
  ctry  <- geo$country
  n_c   <- length(ctry)

  # ---- |-run frame: dates, seed, the source location X ----
  start_date <- as.Date("2025-01-01")
  n_days     <- 300                                   # ~10 months: source rise, export, local waves, peak
  # a large non-EU/EEA source region ("X"): big enough to sustain an epidemic and export travellers
  source <- list(code = "X", name = "Source region X", population = 11e6)

  # ---- |-delay distributions (epidist objects; swap in {epiparameter} via as_epidist()) ----
  # Respiratory-pathogen defaults in the COVID-like range; each is a single, documented assumption.
  delays <- list(
    generation_interval = epidist_gamma("generation interval", mean = 5.2, sd = 1.7),  # Cori-style GI
    incubation          = epidist_lognormal("incubation period", mean = 5.5, sd = 2.3),
    onset_to_death      = epidist_gamma("onset->death", mean = 15.0, sd = 6.6),         # Nishiura CFR delay
    onset_to_report     = epidist_lognormal("onset->report", mean = 4.0, sd = 2.5),     # reporting delay
    onset_to_admission  = epidist_gamma("onset->admission", mean = 7.0, sd = 4.0)       # for the Phase 2 forecast
  )

  # ---- |-reproduction-number trajectories: source X and per-country ----
  # Source: a strong early epidemic, a partial control step, then near-critical -- so exports build,
  # then taper. Countries: a shared template with STAGGERED interventions (spread over ~6 weeks),
  # giving heterogeneous local waves that still share structure. Both are change-point step functions.
  rt_source <- step_schedule(t = c(0, 55, 95), value = c(2.6, 1.4, 0.95))
  set.seed(1)                                          # only to lay out intervention timing; run seed is separate
  intervention_days <- round(stats::runif(n_c, 70, 115))
  rt_country <- make_country_rt(ctry, R_start = 1.9, R_post = 0.8, intervention_days = intervention_days)

  # ---- |-ascertainment: time-varying case detection + a separate, higher death detection ----
  # Case ascertainment ramps UP as testing scales (early cases are missed) -- the classic reason
  # naive CFR is biased early. Deaths are detected far more completely and are treated as constant.
  ascertainment <- list(
    case_rho  = step_schedule(t = c(0, 30, 60, 120), value = c(0.05, 0.15, 0.35, 0.5)),  # rho_case,t
    death_rho = 0.9                                                                       # higher, ~complete
  )

  # ---- |-transmission overdispersion (superspreading): the offspring dispersion k ----
  # Individual-level heterogeneity in onward transmission. `Inf` = a Poisson (homogeneous) process,
  # variance = mean -- the COVID/flu default that keeps runs backward-compatible. A FINITE, small k
  # gives superspreading: most infections transmit to no-one and rare clusters dominate (SARS k ~ 0.16,
  # MERS ~ 0.25, COVID ~ 0.1-0.5). This governs how often an imported chain fizzles out stochastically,
  # which is a defining difference between respiratory pathogens -- see documentation/decisions.md.
  dispersion_k <- Inf                                  # Inf = Poisson; e.g. 0.16 = SARS-like superspreading

  # ---- |-infection fatality: a scalar by default; age structure is an optional refinement ----
  # The DGP is not age-stratified (an extension front), so an age-structured IFR is collapsed to an
  # effective IFR by the population age weights. Provide `ifr_age` (a data.frame age_group/ifr/weight)
  # to use it; otherwise `ifr` is used directly.
  ifr <- 0.006                                         # ~0.6%, COVID-like, all-age effective
  ifr_age <- NULL                                      # e.g. data.frame(age_group=, ifr=, weight=) sums weight to 1

  # ---- |-flights: baseline connectivity ~ population, seasonality, a travel-ban shock ----
  # Baseline daily travellers X->c scale with a power of destination population times country noise;
  # a mild annual seasonality modulates them; a travel ban multiplies all volumes after its date.
  # TRUE volumes drive importation; a noisy OBSERVED copy is what analysts actually see.
  set.seed(2)
  flights <- list(
    baseline_scale   = 40,                             # travellers/day to a 1-million-pop country at connectivity 1
    population_power  = 0.9,                            # volume ~ population^power (larger countries better connected)
    country_noise    = stats::setNames(exp(stats::rnorm(n_c, 0, 0.35)), ctry),  # persistent per-country connectivity
    seasonality_amp  = 0.15,                            # +/-15% annual modulation
    seasonality_peak_day = 210,                         # day of peak travel within the year
    travel_ban_date  = start_date + 100,               # NULL to disable
    travel_ban_factor = 0.15,                           # volumes x this after the ban
    observation_cv   = 0.10                             # lognormal measurement noise on OBSERVED volumes
  )

  # ---- |-per-country surveillance quality (scales detection; drives the importation signal) ----
  # In (0, 1]: 1 = a country that detects (nearly) every imported infection, low = one that misses
  # most. This is what makes the De Salazar importation regression work -- well-surveilled countries
  # reveal the true imports-vs-flights line, and under-detectors fall below it.
  set.seed(3)
  surveillance_quality <- stats::setNames(stats::runif(n_c, 0.3, 1.0), ctry)

  # ---- |-variant: a fitter strain introduced at the source (toggle on/off) ----
  # When enabled, a second strain appears at X on `intro_date` with a multiplicative transmissibility
  # (fitness) advantage; it shares the susceptible pool (complete cross-immunity -- a documented
  # simplification). Imports carry the variant in proportion to the source's variant frequency, so
  # local variant emergence follows from the connectivity, not a separate per-country seeding.
  variant <- list(
    enabled    = TRUE,
    intro_date = start_date + 40,
    fitness    = 0.5,                                  # +50% transmissibility over the wild type
    seed_infections = 5,                               # initial variant infections planted at X
    sequencing_prob = 0.08                             # fraction of reported cases sequenced (typed)
  )

  # ---- |-hospital admissions + healthcare capacity (for the Phase 2 breach check) ----
  # A fraction of infections is admitted (well-ascertained hospital data, dated by onset->admission).
  # Capacity is the ICU/hospital beds available for this pathogen per million population; the Phase 2
  # forecast projects admissions and compares them to capacity_per_million x country population.
  admission_rate       <- 0.02                        # fraction of infections admitted
  capacity_per_million <- 30                           # beds per million (breach threshold)

  list(
    seed         = 2025,                               # THE run seed -- set.seed(seed) inside simulate_pandemic()
    start_date   = start_date,
    n_days       = n_days,
    source       = source,
    countries    = ctry,
    geography    = geo,
    delays       = delays,
    rt_source    = rt_source,
    rt_country   = rt_country,
    dispersion_k = dispersion_k,
    ascertainment = ascertainment,
    ifr          = ifr,
    ifr_age      = ifr_age,
    flights      = flights,
    surveillance_quality = surveillance_quality,
    variant      = variant,
    admission_rate = admission_rate,
    capacity_per_million = capacity_per_million,
    source_seed_infections = 10                        # initial wild-type infections planted at X on day 0
  )
}

# ---- |-restrict a config to a subset of countries (and optionally a shorter run) ----
# Subsets every country-indexed field consistently (geography, per-country R, surveillance quality,
# flight noise), so a smaller experiment -- or a fast test run -- stays internally valid. Editing
# these by hand is the easy way to get a ragged config; this does it in one place.
config_subset <- function(cfg, countries, n_days = NULL) {
  miss <- setdiff(countries, cfg$countries)
  if (length(miss)) stop("config_subset: unknown country/countries: ", paste(miss, collapse = ", "))
  cfg$countries            <- countries
  cfg$geography            <- cfg$geography[cfg$geography$country %in% countries, ]
  cfg$rt_country           <- cfg$rt_country[countries]
  cfg$surveillance_quality <- cfg$surveillance_quality[countries]
  cfg$flights$country_noise <- cfg$flights$country_noise[countries]
  if (!is.null(n_days)) cfg$n_days <- n_days
  cfg
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Validation ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-check a config is internally consistent before a run (fail loud, fail early) ----
# Catches the mistakes that would otherwise surface as a cryptic error deep in the pipeline: missing
# delays, ragged country lists, schedules that don't start at day 0, out-of-range probabilities.
validate_config <- function(cfg) {
  need <- c("seed","start_date","n_days","source","countries","delays","rt_source","rt_country",
            "ascertainment","ifr","flights","surveillance_quality","variant")
  miss <- setdiff(need, names(cfg))
  if (length(miss)) stop("validate_config: missing field(s): ", paste(miss, collapse = ", "))

  if (cfg$n_days < 30) stop("validate_config: n_days too short (< 30) for a recognisable epidemic")

  need_delays <- c("generation_interval","incubation","onset_to_death","onset_to_report")
  miss_d <- setdiff(need_delays, names(cfg$delays))
  if (length(miss_d)) stop("validate_config: missing delay(s): ", paste(miss_d, collapse = ", "))
  for (nm in names(cfg$delays))
    if (!inherits(cfg$delays[[nm]], "epidist"))
      stop(sprintf("validate_config: delay '%s' is not an <epidist> object", nm))

  cc <- cfg$countries
  if (!all(cc %in% names(cfg$surveillance_quality)))   # every country needs a surveillance value
    stop("validate_config: surveillance_quality is missing some countries")
  if (!all(cc %in% names(cfg$rt_country)))
    stop("validate_config: rt_country is missing some countries")

  if (cfg$ascertainment$death_rho <= 0 || cfg$ascertainment$death_rho > 1)
    stop("validate_config: death_rho must be in (0, 1]")
  if (any(cfg$ascertainment$case_rho$value <= 0) || any(cfg$ascertainment$case_rho$value > 1))
    stop("validate_config: case_rho values must be in (0, 1]")
  if (cfg$ifr <= 0 || cfg$ifr >= 1) stop("validate_config: ifr must be in (0, 1)")

  # transmission overdispersion k must be positive (Inf = Poisson); admission rate a probability
  if (!is.null(cfg$dispersion_k) && cfg$dispersion_k <= 0)
    stop("validate_config: dispersion_k must be > 0 (Inf for a Poisson / no-superspreading process)")
  if (!is.null(cfg$admission_rate) && (cfg$admission_rate < 0 || cfg$admission_rate > 1))
    stop("validate_config: admission_rate must be in [0, 1]")

  invisible(TRUE)
}

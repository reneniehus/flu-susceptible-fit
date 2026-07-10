# assemble.R
#
# The top-level pipeline: `simulate_pandemic(config)` runs the whole data-generating process and
# observation model end to end and returns ONE object holding both worlds --
#   $truth    : the latent ground truth (infections, R_t, ascertainment, IFR, imports, flight
#               volumes, variant frequency, source size, prevalence) -- what an analyst never sees;
#   $observed : the degraded surveillance data (reporting triangle, cases by onset & report date,
#               deaths, detected imports + surveillance covariate, sequenced variant counts, observed
#               flight volumes) -- the only thing the analysis toolbox is allowed to use;
#   $config   : the exact config that produced the run (for reproducibility and provenance);
#   $par, $latent : resolved parameters and raw latent stage outputs, for advanced scoring / debugging.
#
# REPRODUCIBILITY. `set.seed(config$seed)` is called once here, so the same config gives byte-identical
# truth and observed data every time. Change the seed for a fresh stochastic realisation; change a
# config field for a different experiment. That is the entire contract.
#
# Requires all of R/ sourced first (see source_all.R).

# ---- |-run the full playground: config -> list(truth, observed, config) ----
simulate_pandemic <- function(cfg = default_config(), quiet = TRUE) {
  validate_config(cfg)
  par <- draw_parameters(cfg)

  set.seed(cfg$seed)                                   # the one place randomness is seeded
  source_x <- simulate_source(cfg, par)
  flights  <- simulate_flights(cfg, par)
  imports  <- simulate_importations(cfg, par, source_x, flights)
  local    <- simulate_local(cfg, par, imports$imports)
  obs      <- observe(cfg, par, source_x, flights, imports, local)

  sim <- structure(
    list(truth = obs$truth, observed = obs$observed, config = cfg, par = par,
         latent = list(source = source_x, flights = flights, imports = imports, local = local)),
    class = "pandemic_sim")

  if (!quiet) print(sim)
  sim
}

# ---- |-a readable one-screen summary of a simulated pandemic ----
print.pandemic_sim <- function(x, ...) {
  cfg <- x$config; tr <- x$truth; ob <- x$observed
  n_ctry <- length(cfg$countries)
  total_inf   <- sum(tr$infections$infections)
  total_cases <- sum(ob$cases_by_report$cases)
  total_death <- sum(ob$deaths$deaths_by_date)
  src_inf     <- sum(tr$source_size$infections)
  peak_src    <- which.max(tr$source_size$infections) - 1
  seeded      <- sum(tapply(tr$imports$imports, tr$imports$country, sum) > 0)

  cat(sprintf("<pandemic_sim>  %s  ->  %d EU/EEA countries   [%d days from %s, seed %d]\n",
              cfg$source$code, n_ctry, cfg$n_days, format(cfg$start_date), cfg$seed))
  cat(sprintf("  variant        : %s\n",
              if (isTRUE(cfg$variant$enabled))
                sprintf("ON  (intro %s, +%.0f%% fitness)", format(cfg$variant$intro_date), 100 * cfg$variant$fitness)
              else "off"))
  cat(sprintf("  source epidemic: %s infections, peak day %d\n", .fmt(src_inf), peak_src))
  cat(sprintf("  countries seeded: %d/%d by importation\n", seeded, n_ctry))
  cat(sprintf("  TRUTH   : %s infections total\n", .fmt(total_inf)))
  cat(sprintf("  OBSERVED: %s reported cases, %s deaths, %d triangle cells\n",
              .fmt(total_cases), .fmt(total_death), nrow(ob$reporting_triangle)))
  cat(sprintf("  effective IFR: %.2f%%   |   case ascertainment: %.0f%% -> %.0f%%\n",
              100 * x$par$ifr_eff, 100 * min(cfg$ascertainment$case_rho$value),
              100 * max(cfg$ascertainment$case_rho$value)))
  invisible(x)
}

# ---- |-compact number formatting for the summary ----
.fmt <- function(x) {
  if (x >= 1e6) sprintf("%.1fM", x / 1e6)
  else if (x >= 1e3) sprintf("%.0fk", x / 1e3)
  else sprintf("%.0f", x)
}

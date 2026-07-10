# analysis_common.R
#
# The seam between the simulation and the analysis toolbox -- and the seam where REAL data enters.
#
# Every analysis tool consumes tidy data frames in a fixed, documented schema, NEVER the raw sim
# object. `as_analysis_input()` pulls those frames out of a simulated pandemic; to analyse real
# surveillance data instead, build a list with the same frames and column names (see
# documentation/real_data.md) and hand it to the very same tools. The tools cannot tell the
# difference -- that is the whole point of routing everything through this schema.
#
# The analyst is assumed to know the DELAY DISTRIBUTIONS from the literature (as they would in a real
# response), so `as_analysis_input()` also carries the config's delays as the analyst's assumed
# delays. Pass your own `epidist` objects to a tool to test sensitivity to delay misspecification --
# a deliberate extension front.
#
# `pp_score()` grades an estimate against the known truth (`sim$truth`). This is what makes the
# playground a testbed rather than a demo: run a tool, score it, see how close it got.
#
# Requires: the engine (R/*) sourced; dplyr available.

suppressMessages({library(dplyr)})

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### The analysis input schema (simulated OR real) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-extract the analyst-facing data (and assumed delays) from a simulated pandemic ----
# Returns the exact schema the tools expect. A real-data user assembles the same list by hand.
as_analysis_input <- function(sim) {
  ob  <- sim$observed
  cfg <- sim$config
  list(
    cases_by_onset     = ob$cases_by_onset,       # location, day, date, cases (right-truncated by onset)
    cases_by_report    = ob$cases_by_report,      # location, day, date, cases (by report date; complete)
    reporting_triangle = ob$reporting_triangle,   # location, onset_day, report_day, cases (nowcast input)
    deaths             = ob$deaths,               # location, day, date, deaths_by_date, deaths_by_onset
    admissions         = ob$admissions,           # location, day, date, admissions (well-ascertained)
    detected_imports   = ob$detected_imports,     # country, day, date, detected_imports, surveillance_quality
    variant_cases      = ob$variant_cases,        # location, day, date, sequenced, variant
    flight_volumes     = ob$flight_volumes,       # country, day, date, volume (observed, noisy)
    delays             = cfg$delays,              # the analyst's assumed (literature) delay distributions
    as_of              = ob$as_of,                # the data cutoff (day index)
    countries          = cfg$countries,
    source_code        = cfg$source$code
  )
}

# ---- |-a location's observed series as a plain day-ordered data frame ----
# Convenience for the tools: pull one location's rows from a schema table, ordered by day.
loc_series <- function(tbl, location) {
  if (is.null(tbl) || !is.data.frame(tbl))
    stop("loc_series: expected a data frame but got ", if (is.null(tbl)) "NULL" else class(tbl)[1],
         " -- a schema table this tool needs is missing from the analysis input. ",
         "See documentation/real_data.md for the required frames.")
  key <- if ("location" %in% names(tbl)) "location" else "country"
  tbl[tbl[[key]] == location, , drop = FALSE] |> (\(d) d[order(d$day), ])()
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Scoring against truth ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-generic point-estimate scorecard: estimate vs truth ----
# Aligns an estimate (optionally with a lower/upper interval) to the truth and returns the standard
# metrics: bias, MAE, RMSE, correlation, and interval coverage. Vectorised over aligned pairs.
pp_score <- function(estimate, truth, lower = NULL, upper = NULL) {
  ok <- is.finite(estimate) & is.finite(truth)
  e <- estimate[ok]; t <- truth[ok]
  out <- list(
    n       = length(e),
    bias    = mean(e - t),
    mae     = mean(abs(e - t)),
    rmse    = sqrt(mean((e - t)^2)),
    rel_bias = mean((e - t) / pmax(abs(t), .Machine$double.eps)),
    cor     = if (length(e) > 2) stats::cor(e, t) else NA_real_
  )
  if (!is.null(lower) && !is.null(upper))
    out$coverage <- mean(truth[ok] >= lower[ok] & truth[ok] <= upper[ok])
  out
}

# ---- |-pretty-print a scorecard as one line ----
pp_score_line <- function(label, s) {
  cov <- if (!is.null(s$coverage)) sprintf("  coverage=%.0f%%", 100 * s$coverage) else ""
  cat(sprintf("%-28s n=%d  bias=%+.3g  MAE=%.3g  RMSE=%.3g  cor=%.2f%s\n",
              label, s$n, s$bias, s$mae, s$rmse, s$cor %||% NA, cov))
  invisible(s)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Truth accessors (what each tool is scored against) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-true reproduction number for a location, by day ----
truth_rt <- function(sim, location) loc_series(sim$truth$Rt, location)

# ---- |-true (eventual detected) cases by onset day for a location ----
truth_cases_by_onset <- function(sim, location) loc_series(sim$truth$cases_by_onset, location)

# ---- |-true infections for a location, by day ----
truth_infections <- function(sim, location) loc_series(sim$truth$infections, location)

# ---- |-the true infection fatality ratio (scalar) ----
truth_ifr <- function(sim) sim$truth$ifr

# ---- |-the true cumulative source-epidemic size, by day ----
truth_source_size <- function(sim) sim$truth$source_size

# ---- |-true imports for a country, by day ----
truth_imports <- function(sim, country) loc_series(sim$truth$imports, country)

# ---- |-true variant frequency for a location, by day ----
truth_variant_freq <- function(sim, location) loc_series(sim$truth$variant_freq, location)

# ---- |-healthcare capacity (admissions threshold) for a country ----
truth_capacity <- function(sim, country) sim$truth$capacity$capacity[sim$truth$capacity$country == country]

# ---- |-null-coalescing helper (mirrors the engine's; harmless if already defined) ----
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

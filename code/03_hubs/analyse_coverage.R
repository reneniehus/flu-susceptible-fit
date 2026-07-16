# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Analyse forecast coverage across weeks, indicators and models ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From the submission-level table (load_forecasts.R) we build the coverage views the
# brief asks for: for each indicator, in which weeks were there forecasts, and by how
# many models (+ the ensemble). Everything keys off `role` from the loader so that
# the ensemble and the reference baseline are always counted apart from real models.

# ---- |-label a round with its winter season ----
# RespiCast rounds run across a winter; anything from August on belongs to the season
# that starts that calendar year (e.g. 2024-10-23 -> "2024/25").
label_season <- function(date) {
  y  <- lubridate::year(date)
  yr <- ifelse(lubridate::month(date) >= 8, y, y - 1L)
  sprintf("%d/%02d", yr, (yr + 1L) %% 100L)
}

# ---- |-the master weekly-coverage table ----------------------------------------
# one row per (indicator, round): how many genuine models, whether the ensemble and
# the baseline were published, and how many countries the ensemble covered.
coverage_weekly <- function(submissions) {
  # the published-ensemble country reach per round (prefer the hub ensemble)
  ens <- submissions %>%
    filter(role == "ensemble") %>%
    group_by(indicator, origin_date) %>%
    summarise(ensemble_locations = max(n_locations), .groups = "drop")

  submissions %>%
    group_by(hub, indicator, origin_date) %>%
    summarise(
      n_models      = sum(role == "model"),         # genuine models feeding the ensemble
      has_ensemble  = any(role == "ensemble"),
      has_baseline  = any(role == "baseline"),
      models_max_locations = safe_max(n_locations[role == "model"]),
      .groups = "drop"
    ) %>%
    left_join(ens, by = c("indicator", "origin_date")) %>%
    mutate(
      season             = label_season(origin_date),
      ensemble_locations = ifelse(is.finite(ensemble_locations), ensemble_locations, NA_integer_),
      # total published forecasters a stakeholder sees that week (models + ensemble)
      n_published        = n_models + as.integer(has_ensemble)
    ) %>%
    arrange(indicator, origin_date)
}

# ---- |-per-model activity ribbon (who was active, when) -------------------------
# one row per (indicator, model): its role, span of rounds and how many it hit.
model_activity <- function(submissions) {
  submissions %>%
    group_by(hub, indicator, model, role) %>%
    summarise(
      n_rounds    = n_distinct(origin_date),
      first_round = min(origin_date),
      last_round  = max(origin_date),
      .groups = "drop"
    ) %>%
    arrange(indicator, role, desc(n_rounds))
}

# the raw presence grid (indicator x model x round) -- feeds the artefact ribbon directly
model_presence <- function(submissions) {
  submissions %>%
    distinct(hub, indicator, model, role, origin_date) %>%
    mutate(season = label_season(origin_date)) %>%
    arrange(indicator, model, origin_date)
}

# ---- |-country coverage --------------------------------------------------------
# union of countries ever forecast, per indicator (splits the per-cell location sets)
country_coverage_overall <- function(submissions) {
  submissions %>%
    filter(role %in% c("model", "ensemble")) %>%
    separate_rows(locations, sep = ",") %>%
    filter(locations != "") %>%
    group_by(indicator) %>%
    summarise(
      n_countries = n_distinct(locations),
      countries   = paste(sort(unique(locations)), collapse = ","),
      .groups = "drop"
    )
}

# ---- |-headline numbers for the artefact ---------------------------------------
coverage_headline <- function(submissions, weekly) {
  per_indicator <- weekly %>%
    group_by(indicator) %>%
    summarise(
      first_round     = min(origin_date),
      last_round      = max(origin_date),
      n_rounds        = n_distinct(origin_date),
      n_seasons       = n_distinct(season),
      peak_models     = max(n_models),
      median_models   = median(n_models),
      rounds_with_ens = sum(has_ensemble),
      .groups = "drop"
    )
  distinct_models <- submissions %>%
    filter(role == "model") %>%
    group_by(indicator) %>%
    summarise(distinct_models = n_distinct(model), .groups = "drop")
  per_indicator %>% left_join(distinct_models, by = "indicator")
}

# per (indicator, season) roll-up for the summary cards
season_summary <- function(submissions, weekly) {
  weekly %>%
    group_by(indicator, season) %>%
    summarise(
      n_rounds      = n_distinct(origin_date),
      first_round   = min(origin_date),
      last_round    = max(origin_date),
      mean_models   = round(mean(n_models), 1),
      peak_models   = max(n_models),
      .groups = "drop"
    ) %>%
    left_join(
      submissions %>%
        filter(role == "model") %>%
        mutate(season = label_season(origin_date)) %>%
        group_by(indicator, season) %>%
        summarise(distinct_models = n_distinct(model), .groups = "drop"),
      by = c("indicator", "season")
    ) %>%
    arrange(indicator, season)
}

# ---- |-run everything and persist ----------------------------------------------
run_coverage_analysis <- function(submissions, params) {
  step("Analysing forecast coverage")
  weekly <- coverage_weekly(submissions)
  out <- list(
    coverage_weekly  = weekly,
    model_activity   = model_activity(submissions),
    model_presence   = model_presence(submissions),
    country_coverage = country_coverage_overall(submissions),
    headline         = coverage_headline(submissions, weekly),
    season_summary   = season_summary(submissions, weekly)
  )

  dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
  # keep the full submissions table too (the reproducible base for everything above)
  write_csv(submissions, file.path(params$output_dir, "hub_submissions.csv"))
  for (nm in names(out)) {
    write_csv(out[[nm]], file.path(params$output_dir, paste0("hub_", nm, ".csv")))
  }
  say(sprintf("wrote submissions + %d coverage tables to output/", length(out)))
  out
}

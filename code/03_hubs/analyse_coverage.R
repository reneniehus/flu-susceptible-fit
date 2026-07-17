# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Analyse forecast coverage across weeks, indicators, models and eras ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From the submission-level table (load_forecasts.R) we build the coverage views: for each indicator,
# in which weeks were there forecasts, by how many models (+ the ensemble), and -- now spanning the
# archived predecessor hubs -- whether the record is CONTINUOUS across the archive -> current handover.
#
# Legacy (Monday) and modern (Wednesday) rounds are snapped to the ISO-week Monday so the two eras land
# on ONE weekly grid; that shared grid is what makes gap-detection (the "uninterrupted since 2021"
# claim) meaningful.

# ---- |-snap any round date to the Monday of its ISO week (the shared weekly key) ----
week_monday <- function(date) lubridate::floor_date(as.Date(date), unit = "week", week_start = 1)

# ---- |-label a round with its winter season (Aug-anchored) ----
label_season <- function(date) {
  y  <- lubridate::year(date)
  yr <- ifelse(lubridate::month(date) >= 8, y, y - 1L)
  sprintf("%d/%02d", yr, (yr + 1L) %% 100L)
}

# ---- |-the master weekly-coverage table (indicator x ISO-week, across eras) ---------------------
# one row per (indicator, week): genuine models, ensemble/baseline presence, ensemble country reach,
# and the era (archive / current) so the handover can be drawn.
coverage_weekly <- function(submissions) {
  s <- submissions %>% mutate(week = week_monday(origin_date))

  ens <- s %>% filter(role == "ensemble") %>%
    group_by(indicator, week) %>% summarise(ensemble_locations = max(n_locations), .groups = "drop")

  s %>%
    group_by(indicator, week) %>%
    summarise(
      n_models     = n_distinct(model[role == "model"]),   # distinct genuine models that week
      has_ensemble = any(role == "ensemble"),
      has_baseline = any(role == "baseline"),
      era          = if (any(era == "current")) "current" else "archive",
      hub          = paste(sort(unique(hub)), collapse = "+"),
      .groups = "drop"
    ) %>%
    left_join(ens, by = c("indicator", "week")) %>%
    mutate(
      season             = label_season(week),
      ensemble_locations = ifelse(is.finite(ensemble_locations), ensemble_locations, NA_integer_),
      n_published        = n_models + as.integer(has_ensemble)
    ) %>%
    arrange(indicator, week)
}

# ---- |-CONTINUITY: is an indicator's weekly record unbroken from first to last round? -----------
# Walks the full Monday grid from the indicator's first to last covered week and finds the gaps.
# Returns one row per indicator with the span, coverage and the longest gap -- the direct test of
# "produced almost uninterrupted since ...".
coverage_continuity <- function(weekly) {
  weekly %>%
    group_by(indicator) %>%
    group_modify(function(d, key) {
      d <- arrange(d, week)
      full  <- seq(min(d$week), max(d$week), by = 7)          # every Monday in the span
      have  <- full %in% d$week
      # run-length encode the missing weeks to size the gaps
      miss  <- !have
      runs  <- rle(miss)
      gap_lengths <- runs$lengths[runs$values]
      tibble(
        first_week    = min(d$week),
        last_week     = max(d$week),
        span_weeks    = length(full),
        covered_weeks = sum(have),
        missing_weeks = sum(miss),
        coverage_pct  = round(100 * sum(have) / length(full), 1),
        n_gaps        = length(gap_lengths),
        longest_gap   = if (length(gap_lengths)) max(gap_lengths) else 0L
      )
    }) %>%
    ungroup()
}

# the actual gap intervals for an indicator (for annotating the timeline)
coverage_gaps <- function(weekly, indicator_name) {
  d <- weekly %>% filter(indicator == indicator_name) %>% arrange(week)
  if (nrow(d) < 2) return(tibble())
  full <- seq(min(d$week), max(d$week), by = 7)
  miss <- !(full %in% d$week)
  r <- rle(miss); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  idx <- which(r$values)
  tibble(indicator = indicator_name,
         gap_start = full[starts[idx]], gap_end = full[ends[idx]],
         gap_weeks = r$lengths[idx]) %>% arrange(desc(gap_weeks))
}

# ---- |-per-model activity + presence (on the ISO-week grid) -------------------------------------
model_activity <- function(submissions) {
  submissions %>% mutate(week = week_monday(origin_date)) %>%
    group_by(hub, era, indicator, model, role) %>%
    summarise(n_rounds = n_distinct(week), first_round = min(week), last_round = max(week), .groups = "drop") %>%
    arrange(indicator, role, desc(n_rounds))
}

model_presence <- function(submissions) {
  submissions %>% mutate(week = week_monday(origin_date)) %>%
    distinct(hub, era, indicator, model, role, week) %>%
    mutate(season = label_season(week)) %>%
    arrange(indicator, model, week)
}

# ---- |-country coverage (union per indicator, across eras) --------------------------------------
country_coverage_overall <- function(submissions) {
  submissions %>%
    filter(role %in% c("model", "ensemble")) %>%
    separate_rows(locations, sep = ",") %>%
    filter(locations != "") %>%
    group_by(indicator) %>%
    summarise(n_countries = n_distinct(locations),
              countries   = paste(sort(unique(locations)), collapse = ","), .groups = "drop")
}

# ---- |-headline numbers per indicator -----------------------------------------------------------
coverage_headline <- function(submissions, weekly, continuity) {
  per <- weekly %>%
    group_by(indicator) %>%
    summarise(first_round = min(week), last_round = max(week),
              n_rounds = n_distinct(week), n_seasons = n_distinct(season),
              peak_models = max(n_models), median_models = median(n_models),
              rounds_with_ens = sum(has_ensemble), .groups = "drop")
  dm <- submissions %>% filter(role == "model") %>%
    group_by(indicator) %>% summarise(distinct_models = n_distinct(model), .groups = "drop")
  per %>% left_join(dm, by = "indicator") %>%
    left_join(select(continuity, indicator, span_weeks, covered_weeks, coverage_pct, n_gaps, longest_gap),
              by = "indicator")
}

# per (indicator, season) roll-up
season_summary <- function(submissions, weekly) {
  weekly %>%
    group_by(indicator, season) %>%
    summarise(n_rounds = n_distinct(week), first_round = min(week), last_round = max(week),
              mean_models = round(mean(n_models), 1), peak_models = max(n_models), .groups = "drop") %>%
    left_join(
      submissions %>% filter(role == "model") %>% mutate(season = label_season(week_monday(origin_date))) %>%
        group_by(indicator, season) %>% summarise(distinct_models = n_distinct(model), .groups = "drop"),
      by = c("indicator", "season")) %>%
    arrange(indicator, season)
}

# ---- |-run everything and persist ---------------------------------------------------------------
run_coverage_analysis <- function(submissions, params) {
  step("Analysing forecast coverage")
  weekly     <- coverage_weekly(submissions)
  continuity <- coverage_continuity(weekly)
  out <- list(
    coverage_weekly  = weekly,
    continuity       = continuity,
    gaps_covid_hosp  = coverage_gaps(weekly, "COVID-19 hospitalisations"),
    model_activity   = model_activity(submissions),
    model_presence   = model_presence(submissions),
    country_coverage = country_coverage_overall(submissions),
    headline         = coverage_headline(submissions, weekly, continuity),
    season_summary   = season_summary(submissions, weekly)
  )

  dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(submissions, file.path(params$output_dir, "hub_submissions.csv"))
  for (nm in names(out)) write_csv(out[[nm]], file.path(params$output_dir, paste0("hub_", nm, ".csv")))
  say(sprintf("wrote submissions + %d coverage tables to output/", length(out)))
  out
}

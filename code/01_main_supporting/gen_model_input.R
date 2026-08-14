# gen_model_input.R
#
# Turns the raw, per-stream data lists (data$epi, data$vax, ...) into a small set of
# tidy, ready-to-plot / ready-to-model tables. Everything is built up from one canonical
# long table; the wide table and the per-country-&-season summary are derived from it.
#
# Outputs of gen_model_input() (the models_in list):
#   - data_timeseries_long : one row per country x season x date x stream x pathogen x
#                            indicator x agegroup; the single source of truth
#   - data_timeseries_wide : same content pivoted so each indicator series is a column
#                            (one row per country x season x date x agegroup)
#   - data_season_summary  : one row per country x season x series with summary / quality
#                            stats (completeness, peak, sums, ...)
#   - contacts             : transformed contact matrices (unchanged from before)
#
# The "+" indicators (ILI+ = ILI consultation rate x pathogen positivity) are built here
# for every pathogen ERVISS reports (Influenza, SARS-CoV-2, RSV), not just influenza.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Small labelling / season helpers ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-season label that a weekly date belongs to (Aug -> Jul window from params) ----
# Compares against the full season-start month-day (not just the month), so it stays consistent
# with add_season_time_columns() for any params$season_start_monthday, not only "-XX-01".
season_of_date = function(date, params=NULL){
  start_monthday  = if (!is.null(params$season_start_monthday)) params$season_start_monthday else "-08-01"
  start_this_year = as.Date(paste0(year(date), start_monthday))
  season_start_year = year(date) - (date < start_this_year)         # before the start date -> previous season
  paste0(season_start_year, "/", season_start_year + 1)
}

# ---- |-human-readable indicator label (pathogen appended where relevant) ----
indicator_label_of = function(indicator, pathogen=NA_character_){
  base = case_when(
    indicator == "ILIconsultationrate" ~ "ILI consultation rate",
    indicator == "ARIconsultationrate" ~ "ARI consultation rate",
    indicator == "ili_plus"            ~ "ILI+",
    indicator == "positivity"          ~ "Positivity",
    indicator == "detections"          ~ "Detections",
    indicator == "tests"               ~ "Tests",
    indicator == "vaccine_coverage"    ~ "Vaccination coverage",
    .default = indicator
  )
  ifelse(is.na(pathogen) | pathogen == "", base, paste0(base, " (", pathogen, ")"))
}

# ---- |-unit attached to each indicator ----
indicator_unit_of = function(indicator){
  case_when(
    indicator %in% c("positivity", "vaccine_coverage")                        ~ "proportion",
    indicator %in% c("detections", "tests")                                   ~ "count",
    indicator %in% c("ILIconsultationrate", "ARIconsultationrate", "ili_plus") ~ "rate",
    .default = NA_character_
  )
}

# ---- |-coerce a stream-specific frame into the canonical long schema ----
# the per-stream extractors below only need to provide country_short, date, indicator, value
# (and optionally agegroup / pathogen / scenario / season); this fills in everything else so
# that every stream lands in exactly the same column layout before they are bound together.
to_canonical_long = function(df, source, stream, params, temporal_resolution="weekly"){
  if (!"agegroup"  %in% names(df)) df$agegroup  = "age_total"
  if (!"pathogen"  %in% names(df)) df$pathogen  = NA_character_
  if (!"scenario"  %in% names(df)) df$scenario  = NA_character_
  if (!"season"    %in% names(df)) df$season    = season_of_date(df$date, params)
  df %>%
    mutate(
      source              = source,
      stream              = stream,
      value               = as.numeric(value),
      indicator_label     = indicator_label_of(indicator, pathogen),
      unit                = indicator_unit_of(indicator),
      temporal_resolution = temporal_resolution,
      observed            = !is.na(value)
    ) %>%
    add_season_time_columns(params=params) %>%
    select(country_short, season, date,
           season_start_year, season_start_date, season_day, season_week, iso_week,
           source, stream, pathogen, scenario, indicator, indicator_label, agegroup,
           value, unit, temporal_resolution, observed)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Per-stream extractors (raw data lists -> canonical long pieces) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-ERVISS ILI / ARI consultation rates (age-specific) ----
extract_iliari = function(data, params){
  data$epi$erviss_ili_ari %>%                                  # country_short, date, target, agegroup, value
    transmute(country_short, date=as.Date(date), agegroup, indicator=target, value) %>%
    to_canonical_long(source="ERVISS", stream="ili_ari", params=params)
}

# ---- |-ERVISS typing: tests, detections and (re-computed) positivity per pathogen ----
# tests / detections come straight from the file; positivity is re-derived as detections/tests
# (robust to gaps in the published positivity column). Typing is only reported at age_total.
# The pathogen-level total row (the one carrying both detections AND the shared tests
# denominator) is the row where pathogentype == pathogen: "Influenza"/total, "RSV"/RSV,
# "SARS-CoV-2"/SARS-CoV-2 -- the finer pathogentype/subtype rows only carry detections.
extract_typing = function(epi_key, stream, data, params){
  check_required_columns(data$epi[[epi_key]],                        # typing needs the pathogen breakdown
                         c("pathogen", "pathogentype", "pathogensubtype", "age", "indicator", "value"),
                         paste0("data$epi$", epi_key))
  wide = data$epi[[epi_key]] %>%
    filter(pathogentype == pathogen, age == "age_total",
           indicator %in% c("detections", "tests")) %>%
    transmute(country_short, date=as.Date(date), pathogen, indicator, value=as.numeric(value)) %>%
    # one value per country x date x pathogen x indicator (guards against any duplicate rows)
    summarise(value=sum(value, na.rm=FALSE), .by=c(country_short, date, pathogen, indicator)) %>%
    pivot_wider(names_from=indicator, values_from=value)
  wide %>%
    mutate(positivity = replace_inf(detections / tests, NA)) %>%   # own positivity = detections / tests
    pivot_longer(c(detections, tests, positivity), names_to="indicator", values_to="value") %>%
    to_canonical_long(source="ERVISS", stream=stream, params=params)
}

# ---- |-"+" indicators: syndromic rate x pathogen positivity (ILI+ for every pathogen) ----
# positivity is age_total only, so it is broadcast across the age-specific ILI rate. Produces
# one ili_plus series per pathogen present in the typing stream (Influenza / SARS-CoV-2 / RSV).
build_ili_plus = function(iliari_long, typing_long, stream, params){
  ili = iliari_long %>%
    filter(indicator == "ILIconsultationrate") %>%
    select(country_short, date, agegroup, ili=value)
  positivity = typing_long %>%
    filter(indicator == "positivity") %>%
    select(country_short, date, pathogen, positivity=value)
  ili %>%
    inner_join(positivity, by=c("country_short", "date"),   # broadcast pathogen positivity across ages
               relationship="many-to-many") %>%
    transmute(country_short, date, agegroup, pathogen,
              indicator="ili_plus", value=ili * positivity) %>%
    to_canonical_long(source="ERVISS", stream=stream, params=params)
}

# ---- |-RespiCompass ILI+ (externally provided, influenza only) ----
extract_respicompass_iliplus = function(data, params){
  data$epi$respicompass_iliplus %>%                           # country_short, date, target, agegroup, value
    transmute(country_short, date=as.Date(date), agegroup,
              pathogen="Influenza", indicator="ili_plus", value) %>%
    to_canonical_long(source="RespiCompass", stream="ili_plus_respicompass", params=params)
}

# ---- |-Vaccination coverage: observed history (65+) and forward scenarios (seasonal) ----
extract_vax = function(data, params){
  pieces = list()
  if (!is.null(data$vax$data_vax_history)) {
    pieces$history = data$vax$data_vax_history %>%
      transmute(country_short=iso2_code, season,
                date=as.Date(NA), agegroup=target_group,
                scenario="observed_history",
                indicator="vaccine_coverage", value=as.numeric(vaccine_coverage)) %>%
      to_canonical_long(source="RespiCompass", stream="vaccination_history_65plus",
                        params=params, temporal_resolution="seasonal")
  }
  if (!is.null(data$vax$data_vax)) {
    # the scenarios target the season of the RespiCompass round they were fetched from
    # ("2024-2025_round_1" -> "2024/2025") -- parsed from the round the data actually came from,
    # not stamped from a separate latest-season knob (which historically drifted out of sync
    # and mislabelled this stream by one season)
    round_season = sub("^(\\d{4})-(\\d{4}).*$", "\\1/\\2", params$respicompass_round)
    pieces$scenarios = data$vax$data_vax %>%
      pivot_longer(any_of(c("higher_vax_coverage", "lower_vax_coverage", "no_vaccination")),
                   names_to="scenario", values_to="vaccine_coverage") %>%
      transmute(country_short=iso2_code,
                season=round_season,
                date=as.Date(NA), agegroup=target_group, scenario,
                indicator="vaccine_coverage", value=as.numeric(vaccine_coverage)) %>%
      to_canonical_long(source="RespiCompass", stream="vaccination_scenario",
                        params=params, temporal_resolution="seasonal")
  }
  bind_rows(pieces)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Output builders: long, wide, season summary ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-canonical long table: every indicator stacked in one tidy frame ----
make_data_timeseries_long = function(data=NULL, params=NULL){
  iliari        = extract_iliari(data, params)
  typing_sent   = extract_typing("erviss_typing_sentinel",    "typing_sentinel",    data, params)
  typing_nonsen = extract_typing("erviss_typing_nonsentinel", "typing_nonsentinel", data, params)
  bind_rows(
    iliari,
    typing_sent,
    typing_nonsen,
    build_ili_plus(iliari, typing_sent,   "ili_plus_sentinel",    params),  # ILI+ per pathogen (sentinel)
    build_ili_plus(iliari, typing_nonsen, "ili_plus_nonsentinel", params),  # ILI+ per pathogen (non-sentinel)
    extract_respicompass_iliplus(data, params),                             # ILI+ flu (RespiCompass)
    extract_vax(data, params)
  ) %>%
    arrange(country_short, season, source, stream, pathogen, indicator, scenario, agegroup, date)
}

# ---- |-wide table: one column per indicator series (weekly streams only) ----
# series_id glues stream + indicator + pathogen into a single column name, so the wide table
# is keyed purely by country x season x date x agegroup.
make_data_timeseries_wide = function(data_timeseries_long){
  data_timeseries_long %>%
    filter(temporal_resolution == "weekly") %>%
    mutate(series_id = paste(stream, indicator, sep="__"),
           series_id = ifelse(is.na(pathogen), series_id, paste(series_id, pathogen, sep="__")),
           series_id = make.names(series_id)) %>%
    select(country_short, season, date, season_week, agegroup, series_id, value) %>%
    pivot_wider(names_from=series_id, values_from=value,
                values_fn=function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm=TRUE)) %>%
    arrange(country_short, season, agegroup, date)
}

# ---- |-summary / quality stats for one group of long rows ----
summarise_timeseries_group = function(df){
  # min/max guarded so all-NA date groups (e.g. seasonal vaccination) return NA, not a warning
  first_obs_date = function(date, observed){ d = date[observed & !is.na(date)]; if (length(d)==0) as.Date(NA) else min(d) }
  last_obs_date  = function(date, observed){ d = date[observed & !is.na(date)]; if (length(d)==0) as.Date(NA) else max(d) }
  df %>% summarise(
    n_rows            = n(),
    n_observed        = sum(observed, na.rm=TRUE),
    n_weeks_observed  = n_distinct(date[observed & !is.na(date)]),
    observed_fraction = mean(observed, na.rm=TRUE),
    sum_value         = ifelse(all(is.na(value)), NA_real_, sum(value, na.rm=TRUE)),   # NA-guarded like mean/max
    mean_value        = ifelse(all(is.na(value)), NA_real_, mean(value, na.rm=TRUE)),
    max_value         = ifelse(all(is.na(value)), NA_real_, max(value, na.rm=TRUE)),
    peak_date         = if (all(is.na(value)) || all(is.na(date))) as.Date(NA) else date[which.max(replace_na(value, -Inf))],
    first_date        = first_obs_date(date, observed),
    last_date         = last_obs_date(date, observed),
    .groups="drop"
  )
}

# ---- |-per-country-&-season summary, both age-specific and pooled over age ----
# completeness = observed weeks / number of ISO weeks in that season (weekly streams only).
make_data_season_summary = function(data_timeseries_long, params=NULL){
  # expected number of weekly (Wednesday) data points per season
  weeks_in_season = data_timeseries_long %>%
    distinct(season) %>%
    mutate(
      .season_start_year   = season_start_year_from_label(season),
      .start = ymd(paste0(.season_start_year,     params$season_start_monthday)),
      .end   = ymd(paste0(.season_start_year + 1, params$season_end_monthday)),
      # "%u" = ISO weekday number (Wed = 3): locale-proof, unlike weekdays() == "Wednesday",
      # which returns 0 matches (-> completeness = Inf) under any non-English LC_TIME
      weeks_in_season = map2_int(.start, .end, ~ sum(format(seq(.x, .y, by="day"), "%u") == "3"))
    ) %>%
    select(season, weeks_in_season)

  group_core = c("country_short", "season", "source", "stream", "pathogen",
                 "indicator", "indicator_label", "unit", "temporal_resolution")
  per_age = data_timeseries_long %>%
    group_by(across(all_of(c(group_core, "scenario", "agegroup")))) %>%
    summarise_timeseries_group() %>%
    mutate(summary_level="agegroup")
  pooled_age = data_timeseries_long %>%
    group_by(across(all_of(c(group_core, "scenario")))) %>%
    summarise_timeseries_group() %>%
    mutate(agegroup="all_agegroups", summary_level="all_agegroups")

  # Two complementary completeness measures (a week is "observed" if its value is non-missing;
  # a reported 0 counts as observed, and positivity needs both detections and tests present):
  #   completeness        = observed weeks / ISO weeks in the FULL Aug-Jul season window.
  #                         Penalises the off-season, so winter-only streams cap well below 1.
  #   completeness_active = observed weeks / weeks between the first and last observed week.
  #                         Ignores the off-season; measures gaps WITHIN the active reporting span.
  bind_rows(per_age, pooled_age) %>%
    left_join(weeks_in_season, by="season") %>%
    mutate(
      weeks_in_span       = ifelse(temporal_resolution == "weekly" & !is.na(first_date) & !is.na(last_date),
                                   as.integer(last_date - first_date) %/% 7L + 1L, NA_integer_),
      completeness        = ifelse(temporal_resolution == "weekly",
                                   n_weeks_observed / weeks_in_season, NA_real_),
      completeness_active = ifelse(!is.na(weeks_in_span) & weeks_in_span > 0,
                                   pmin(1, n_weeks_observed / weeks_in_span), NA_real_)
    ) %>%
    select(country_short, season, source, stream, pathogen, indicator, indicator_label,
           scenario, agegroup, summary_level, unit, temporal_resolution,
           weeks_in_season, weeks_in_span, n_weeks_observed,
           completeness, completeness_active, everything()) %>%
    arrange(country_short, season, source, stream, pathogen, indicator, summary_level, agegroup)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### gen_model_input(): assemble the models_in list ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gen_model_input = function( params=NULL , data=NULL ){
  t1 <- Sys.time()

  models_in = list(
    time_of_execution = now(),    # time-stamp
    duration = NULL               # execution duration
  )

  ## ---- |-Canonical long time-series (single source of truth) ----
  models_in$data_timeseries_long = make_data_timeseries_long(data=data, params=params)
  ## ---- |-Wide version (one column per indicator series) ----
  models_in$data_timeseries_wide = make_data_timeseries_wide(models_in$data_timeseries_long)
  ## ---- |-Per-country-&-season summary / quality stats ----
  models_in$data_season_summary  = make_data_season_summary(models_in$data_timeseries_long, params=params)
  ## ---- |-Contacts ----
  models_in$contacts = transform_contacts(data, params)   # transform the contact matrices for model requirements

  #### output
  t2 <- Sys.time()
  models_in$duration = get_in_hms(t2, t1)
  return(models_in)
}

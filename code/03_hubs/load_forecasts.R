# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Load the forecast submissions from the forecasting hubs (current + archived) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Five hubs, two file formats (see config.R $hubs):
#   MODERN (hubverse): RespiCast-Covid19, RespiCast-SyndromicIndicators, flu/ari archives
#     model-output/<team>/<origin_date>-<team>.csv, with an `origin_date` and a `target` column.
#   LEGACY (old EU COVID hub): covid19-forecast-hub-europe_archive
#     data-processed/<team>/<forecast_date>-<team>.csv, with `forecast_date` and compound targets
#     like "2 wk ahead inc hosp" -- the indicator is the last token (case / hosp / death).
#
# We don't care about the forecast VALUES here -- only coverage: who forecast which indicator, for
# which weekly round, for how many countries/horizons. The unit of the returned table is one
# (hub, model, origin_date, indicator) cell, with `era` (archive / current) carried through so the
# handover between an archived hub and its live successor can be told.
#
# Two quirks handled below, both found by inspecting the raw files:
#   1. the modern CSV COLUMN ORDER is not constant across files -> always select by NAME.
#   2. a few submissions are header-only (no data rows) -> skipped, and counted.

# ---- |-classify a model-output folder: model / ensemble / baseline ----
model_role <- function(model, params) {
  if (model %in% params$ensemble_models) return("ensemble")
  if (model %in% params$baseline_models) return("baseline")
  "model"
}

# ---- |-read ONE modern (hubverse) forecast file down to its coverage index ----
read_forecast_index <- function(file, params) {
  dt <- tryCatch(
    data.table::fread(file, select = c("origin_date", "target", "location", "horizon"),
                      showProgress = FALSE, colClasses = list(character = "origin_date")),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  as_tibble(dt) %>%
    mutate(indicator = recode(as.character(target), !!!params$target_labels, .default = as.character(target))) %>%
    group_by(origin_date, indicator) %>%
    summarise(
      n_locations = n_distinct(location),
      n_horizons  = n_distinct(horizon),
      n_rows      = n(),
      locations   = paste(sort(unique(as.character(location))), collapse = ","),
      .groups = "drop"
    )
}

# ---- |-read ONE legacy (old EU COVID hub) file down to its coverage index ----
# targets look like "N wk ahead {inc|cum} {case|hosp|death}"; pull the indicator token and horizon.
read_legacy_index <- function(file, params) {
  dt <- tryCatch(
    data.table::fread(file, select = c("forecast_date", "target", "location"),
                      showProgress = FALSE, colClasses = list(character = "forecast_date")),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  as_tibble(dt) %>%
    mutate(
      token   = str_extract(as.character(target), "(case|hosp|death)$"),
      horizon = as.integer(str_extract(as.character(target), "^-?\\d+"))
    ) %>%
    filter(!is.na(token)) %>%                                   # drop any non-case/hosp/death target
    mutate(indicator   = recode(token, !!!params$target_labels),
           origin_date = as.character(forecast_date)) %>%
    group_by(origin_date, indicator) %>%
    summarise(
      n_locations = n_distinct(location),
      n_horizons  = n_distinct(horizon),
      n_rows      = n(),
      locations   = paste(sort(unique(as.character(location))), collapse = ","),
      .groups = "drop"
    )
}

# ---- |-load every submission from one hub (dispatches on format) ----
load_one_hub <- function(hub_row, params) {
  legacy   <- identical(hub_row$format, "legacy")
  sub_dir  <- if (legacy) "data-processed" else "model-output"
  hub_path <- file.path(params$hubs_dir, hub_row$folder, sub_dir)
  if (!dir.exists(hub_path))
    stop(sprintf("Hub folder not found: %s\n  clone the hubs into %s (see README).", hub_path, params$hubs_dir))

  files <- list.files(hub_path, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  say(sprintf("%s (%s, %s): %d forecast files", hub_row$name, hub_row$format, hub_row$era, length(files)))

  n_empty <- 0L
  read_idx <- if (legacy) read_legacy_index else read_forecast_index
  rows <- map(files, function(f) {
    model <- basename(dirname(f))
    role  <- model_role(model, params)          # scalar, before the mutate (mutate re-reads `model` as a column)
    idx   <- read_idx(f, params)
    if (is.null(idx)) { n_empty <<- n_empty + 1L; return(NULL) }
    idx %>% mutate(hub = hub_row$name, era = hub_row$era, model = model, role = role, file = basename(f))
  })

  out <- bind_rows(rows)
  if (n_empty > 0) say(sprintf("  (skipped %d header-only / unreadable file(s))", n_empty))
  out
}

# ---- |-mother loader across all hubs ----
load_hub_forecasts <- function(params) {
  step("Loading forecast submissions from the forecasting hubs (current + archived)")

  submissions <- params$hubs %>%
    split(seq_len(nrow(.))) %>%
    map(~ load_one_hub(.x, params)) %>%
    bind_rows() %>%
    mutate(origin_date = as.Date(origin_date)) %>%
    select(hub, era, indicator, origin_date, model, role, n_locations, n_horizons, n_rows, locations, file) %>%
    arrange(indicator, origin_date, model)

  say(sprintf("total submission cells: %d | models: %d | rounds: %d | indicators: %d | span %s -> %s",
              nrow(submissions),
              n_distinct(submissions$model),
              n_distinct(submissions$origin_date),
              n_distinct(submissions$indicator),
              min(submissions$origin_date), max(submissions$origin_date)))
  submissions
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Load the forecast submissions from the two RespiCast hubs ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The two hubs follow the hubverse layout: model-output/<team-model>/<origin_date>-<team-model>.csv,
# each file holding one team's probabilistic forecast (quantiles) for a given round.
# We don't care about the forecast VALUES here -- only about coverage: who forecast
# what indicator, for which round (week), and for how many countries/horizons.
#
# Two quirks handled below, both found by inspecting the raw files:
#   1. the CSV COLUMN ORDER is not constant across files -> always select by NAME.
#   2. a few submissions are header-only (no data rows) -> skipped, and counted.
#
# The unit of the returned table is one (hub, model, origin_date, indicator) cell.

# ---- |-classify a model-output folder: model / ensemble / baseline ----
model_role <- function(model, params) {
  if (model %in% params$ensemble_models) return("ensemble")
  if (model %in% params$baseline_models) return("baseline")
  "model"
}

# ---- |-read ONE forecast file down to its coverage index ----
# returns one row per (origin_date, indicator) with country / horizon counts,
# or NULL for an unreadable or header-only file.
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
      locations   = paste(sort(unique(as.character(location))), collapse = ","),  # the country set, for union/coverage
      .groups = "drop"
    )
}

# ---- |-load every submission from one hub ----
load_one_hub <- function(hub_row, params) {
  hub_path <- file.path(params$hubs_dir, hub_row$folder, "model-output")
  if (!dir.exists(hub_path))
    stop(sprintf("Hub folder not found: %s\n  clone the hubs into %s (see README).", hub_path, params$hubs_dir))

  files <- list.files(hub_path, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  say(sprintf("%s: %d forecast files under model-output/", hub_row$name, length(files)))

  n_empty <- 0L
  rows <- map(files, function(f) {
    model <- basename(dirname(f))
    role  <- model_role(model, params)      # scalar, computed before the mutate (mutate would re-read `model` as a column)
    idx   <- read_forecast_index(f, params)
    if (is.null(idx)) { n_empty <<- n_empty + 1L; return(NULL) }
    idx %>% mutate(hub = hub_row$name, model = model, role = role, file = basename(f))
  })

  out <- bind_rows(rows)
  if (n_empty > 0) say(sprintf("  (skipped %d header-only / unreadable file(s))", n_empty))
  out
}

# ---- |-mother loader across both hubs ----
load_hub_forecasts <- function(params) {
  step("Loading forecast submissions from the two RespiCast hubs")

  submissions <- params$hubs %>%
    split(seq_len(nrow(.))) %>%
    map(~ load_one_hub(.x, params)) %>%
    bind_rows() %>%
    mutate(origin_date = as.Date(origin_date)) %>%
    # tidy final column order
    select(hub, indicator, origin_date, model, role, n_locations, n_horizons, n_rows, locations, file) %>%
    arrange(hub, indicator, origin_date, model)

  say(sprintf("total submission cells: %d | models: %d | rounds: %d | indicators: %d",
              nrow(submissions),
              n_distinct(submissions$model),
              n_distinct(submissions$origin_date),
              n_distinct(submissions$indicator)))
  submissions
}

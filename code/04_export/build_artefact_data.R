# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Assemble the artefact JSON ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The HTML dashboard renders entirely client-side from ONE embedded JSON blob (the
# artefact sandbox blocks network fetches, so nothing is loaded at runtime). This
# script gathers the survey + coverage results into that blob and writes it to
# output/artefact_data.json, which code/05_artefact/build_pages.R injects into the
# page template. Keep it a thin assembler: all real computation lives in the
# analyse_*.R modules; here we only shape and label for display.
#
# The coverage half is reshaped from tidy records into arrays aligned to ONE master
# week axis -- the shape a heatmap / line chart wants -- plus a compact per-model
# presence index for the participation ribbon.

# a data frame -> list-of-records (the shape JS consumes for tables / small sets).
# Ordered factors must go out as their LABELS, not their integer codes -- otherwise
# jsonlite serialises e.g. the staff band as 1/2/4 and the front-end can't match it.
records <- function(df) {
  df <- as.data.frame(df)
  is_fac <- vapply(df, is.factor, logical(1))
  df[is_fac] <- lapply(df[is_fac], as.character)
  purrr::transpose(as.list(df))
}

# ---- |-coverage reshaped for the charts ----------------------------------------
# weeks : the sorted union of every round date (the x axis)
# series: per indicator, vectors aligned to `weeks` (n_models, ensemble flag/reach, season)
# ribbon: per (model, indicator), the week indices it was present (0-based, for JS)
build_coverage_web <- function(coverage_out, params) {
  weekly <- coverage_out$coverage_weekly
  weeks  <- sort(unique(as.character(weekly$origin_date)))
  wk_idx <- setNames(seq_along(weeks) - 1L, weeks)          # 0-based index for JS

  indicators <- sort(unique(weekly$indicator))
  series <- lapply(indicators, function(ind) {
    d <- weekly %>% filter(indicator == ind)
    align <- function(col, fill = NA) {
      v <- rep(fill, length(weeks)); v[wk_idx[as.character(d$origin_date)] + 1L] <- d[[col]]; v
    }
    list(
      n_models           = align("n_models", 0L),
      has_ensemble       = as.integer(align("has_ensemble", FALSE)),
      ensemble_locations = align("ensemble_locations"),
      season             = align("season")
    )
  })
  names(series) <- indicators

  ribbon <- coverage_out$model_presence %>%
    group_by(indicator, model, role) %>%
    summarise(present = list(unname(wk_idx[as.character(origin_date)])), .groups = "drop") %>%
    arrange(indicator, role, model)

  list(
    weeks            = weeks,
    indicators       = indicators,
    series           = series,
    ribbon           = records(ribbon %>% mutate(present = present)),
    country_coverage = records(coverage_out$country_coverage),
    headline         = records(coverage_out$headline %>%
                                 mutate(first_round = as.character(first_round),
                                        last_round  = as.character(last_round))),
    season_summary   = records(coverage_out$season_summary %>%
                                 mutate(first_round = as.character(first_round),
                                        last_round  = as.character(last_round)))
  )
}

# ---- |-mother assembler --------------------------------------------------------
build_artefact_data <- function(survey, survey_out, coverage_out, params, themes = NULL) {
  step("Building artefact JSON")

  artefact <- list(
    meta = list(
      generated_from = "code/00_main.R",
      n_respondents  = nrow(survey$responses),
      hubs           = records(params$hubs %>% select(name, folder, indicators))
    ),
    survey = list(
      codebook              = records(survey$codebook %>% select(variable, question, type, label)),
      engagement            = records(survey_out$engagement),
      engagement_dist       = records(survey_out$engagement_dist),
      engagement_capacity   = records(survey_out$engagement_capacity),
      staff                 = records(survey_out$staff),
      decisions_dist        = records(survey_out$decisions_dist %>% mutate(level = as.character(level))),
      decisions_headline    = records(survey_out$decisions_headline),
      value_choice          = records(survey_out$value_choice),
      value_choice_capacity = records(survey_out$value_choice_capacity),
      integration           = records(survey_out$integration %>% mutate(level = as.character(level))),
      rankings              = records(survey_out$rankings),
      free_text             = records(survey_out$free_text)
    ),
    coverage = build_coverage_web(coverage_out, params),
    # qualitative theme synthesis (from the verification/coding workflow); NULL until injected
    themes = themes
  )

  dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(params$output_dir, "artefact_data.json")
  write_json(artefact, path, auto_unbox = TRUE, pretty = FALSE, na = "null")
  say(sprintf("wrote %s (%.0f KB)", path, file.info(path)$size / 1024))
  invisible(artefact)
}

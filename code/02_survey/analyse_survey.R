# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Analyse the NFP survey ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Turns the tidy $responses table into the small set of summary tables the artefact
# and the findings doc consume. Every function takes the survey list from
# load_survey() and returns a plain data frame; run_survey_analysis() calls them
# all and also writes each to output/ as CSV.
#
# Framing: the brief is the *value (and potential value) of forecasting/nowcasting*
# to external stakeholders. So the analysis keeps circling one axis -- in-house
# modelling capacity (Q6) -- because that is what decides whether ECDC forecasting
# is a convenience or the ONLY forecasting a country has.

# the five ECDC outputs, classified so we can compare forecasting vs scenario vs guidance
ENGAGEMENT_CLASS <- tribble(
  ~variable,              ~output,                          ~class,
  "eng_covid_guidance",   "COVID-19 risk assessments",      "Guidance",
  "eng_covid_forecast",   "COVID-19 Forecast Hub",          "Short-term forecasting",
  "eng_respicast",        "RespiCast",                      "Short-term forecasting",
  "eng_covid_scenario",   "COVID-19 Scenario Hub",          "Scenario modelling",
  "eng_respicompass",     "RespiCompass",                   "Scenario modelling"
)

# the five ranked communication channels (Q11/Q12), canonical spellings
COMMS_CHANNELS <- c(
  "Dedicated newsletter / mailing list",
  "Dedicated website",
  "Online meetings for data analysts and modellers",
  "Online meetings for policy decision-makers",
  "Other"
)

# ---- |-Q5: AWARENESS of each ECDC output ----------------------------------------
# NB the function/table names say "engagement" for historical reasons, but Q5 measured
# AWARENESS ("How aware are you...", 0 = not aware, 5 = fully aware) -- NOT use. Every
# human-facing label reads "awareness". Per output: response count, how many gave 0
# (= not aware at all), mean & median on the 0-5 scale, and the share "high" (>=4, fully aware).
engagement_summary <- function(survey) {
  survey$responses %>%
    select(all_of(ENGAGEMENT_CLASS$variable)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "score") %>%
    filter(!is.na(score)) %>%
    group_by(variable) %>%
    summarise(
      n          = n(),
      n_not_aware = sum(score == 0),
      pct_not_aware = round(100 * mean(score == 0)),
      mean_score = round(mean(score), 2),
      median_score = median(score),
      pct_high   = round(100 * mean(score >= 4)),
      .groups = "drop"
    ) %>%
    left_join(ENGAGEMENT_CLASS, by = "variable") %>%
    select(output, class, n, n_not_aware, pct_not_aware, mean_score, median_score, pct_high) %>%
    arrange(desc(mean_score))
}

# full 0-5 distribution per output (for the artefact's small-multiple / heat strip)
engagement_distribution <- function(survey) {
  survey$responses %>%
    select(all_of(ENGAGEMENT_CLASS$variable)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "score") %>%
    filter(!is.na(score)) %>%
    count(variable, score, name = "n") %>%
    left_join(ENGAGEMENT_CLASS, by = "variable") %>%
    select(output, class, score, n) %>%
    arrange(output, score)
}

# ---- |-The capacity axis: engagement & value by in-house modelling staff ---------
# The value story. Mean engagement with each output, split by Q6 staff band, so we
# can see whether the countries that CANNOT model in-house lean on ECDC forecasts
# more or less than the countries that can.
engagement_by_capacity <- function(survey) {
  survey$responses %>%
    select(staff, all_of(ENGAGEMENT_CLASS$variable)) %>%
    pivot_longer(-staff, names_to = "variable", values_to = "score") %>%
    filter(!is.na(score), !is.na(staff)) %>%
    left_join(ENGAGEMENT_CLASS, by = "variable") %>%
    group_by(class, staff) %>%
    summarise(mean_score = round(mean(score), 2), n = n(), .groups = "drop") %>%
    arrange(class, staff)
}

staff_summary <- function(survey) {
  survey$responses %>%
    filter(!is.na(staff)) %>%
    count(staff, name = "n") %>%
    mutate(pct = round(100 * n / sum(n)))
}

# ---- |-Q7: which public-health decisions modelling would inform ------------------
# Long table (action x Likert level x n/pct) feeding a diverging stacked bar, plus
# the headline "share Likely / Very likely" per action.
DECISION_VARS <- tribble(
  ~variable,             ~action,
  "dec_surveillance",    "Inform surveillance activities",
  "dec_vaccination",     "Planning vaccination campaigns",
  "dec_healthcare",      "Healthcare capacity planning",
  "dec_countermeasures", "Procuring medical countermeasures",
  "dec_other",           "Other activities"
)

decisions_distribution <- function(survey, params) {
  survey$responses %>%
    select(all_of(DECISION_VARS$variable)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "level") %>%
    filter(!is.na(level)) %>%
    left_join(DECISION_VARS, by = "variable") %>%
    mutate(level = factor(level, levels = params$levels_q7, ordered = TRUE)) %>%
    count(action, level, name = "n") %>%
    group_by(action) %>%
    mutate(pct = round(100 * n / sum(n))) %>%
    ungroup() %>%
    arrange(action, level)
}

decisions_headline <- function(survey) {
  survey$responses %>%
    select(all_of(DECISION_VARS$variable)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "level") %>%
    filter(!is.na(level)) %>%
    left_join(DECISION_VARS, by = "variable") %>%
    group_by(action) %>%
    summarise(
      n            = n(),
      pct_likely   = round(100 * mean(level %in% c("Likely", "Very likely"))),
      pct_unlikely = round(100 * mean(level %in% c("Unlikely", "Very unlikely"))),
      .groups = "drop"
    ) %>%
    arrange(desc(pct_likely))
}

# ---- |-Q8: RespiCast (forecasts) vs RespiCompass (scenarios) ---------------------
value_choice_summary <- function(survey) {
  survey$responses %>%
    filter(!is.na(value_choice)) %>%
    count(value_choice, name = "n") %>%
    mutate(pct = round(100 * n / sum(n))) %>%
    arrange(desc(n))
}

# ...and the same split by capacity band (does forecasting appeal more where there is no in-house team?)
value_choice_by_capacity <- function(survey) {
  survey$responses %>%
    filter(!is.na(value_choice), !is.na(staff)) %>%
    count(staff, value_choice, name = "n")
}

# ---- |-Q10: is there a mechanism to integrate modelling into decisions? ----------
integration_summary <- function(survey, params) {
  survey$responses %>%
    filter(!is.na(integration)) %>%
    mutate(level = factor(as.character(integration), levels = params$levels_q10, ordered = TRUE)) %>%
    count(level, name = "n") %>%
    mutate(pct = round(100 * n / sum(n)))
}

# ---- |-Q11/Q12: preferred communication channels --------------------------------
# The stored value is a single string with the channels in preference order,
# separated by "; ". Split it, and the position in the list IS the rank (1 = best).
parse_ranking <- function(ranking_string) {
  parts <- str_squish(str_split(ranking_string, ";")[[1]])
  parts <- parts[parts != ""]
  tibble(channel = parts, rank = seq_along(parts))
}

rankings_summary <- function(survey) {
  one_question <- function(col, mode) {
    survey$responses %>%
      filter(!is.na(.data[[col]])) %>%
      transmute(respondent_id, ranking = .data[[col]]) %>%
      rowwise() %>%
      reframe(respondent_id, parse_ranking(ranking)) %>%
      # keep only the five canonical channels (defensive against spelling drift)
      filter(channel %in% COMMS_CHANNELS) %>%
      group_by(channel) %>%
      summarise(
        n           = n(),
        mean_rank   = round(mean(rank), 2),   # lower = more preferred
        pct_first   = round(100 * mean(rank == 1)),
        pct_top2    = round(100 * mean(rank <= 2)),
        .groups = "drop"
      ) %>%
      mutate(mode = mode) %>%
      arrange(mean_rank)
  }
  bind_rows(
    one_question("comms_forecast_rank", "Short-term forecasts (RespiCast, weekly)"),
    one_question("comms_scenario_rank", "Scenario insights (RespiCompass, 1-2x/yr)")
  )
}

# ---- |-Free text: every non-empty open response, tagged by question --------------
free_text_table <- function(survey) {
  text_vars <- c(reflections = "Q13", value_choice_text = "Q8b",
                 dec_other_text = "Q7b", comms_forecast_text = "Q11b",
                 comms_scenario_text = "Q12b")
  survey$responses %>%
    select(respondent_id, all_of(names(text_vars))) %>%
    pivot_longer(-respondent_id, names_to = "variable", values_to = "text") %>%
    filter(!is.na(text), str_length(text) > 2) %>%
    mutate(question = text_vars[variable]) %>%
    select(respondent_id, question, variable, text)
}

# ---- |-Run everything and persist -----------------------------------------------
run_survey_analysis <- function(survey, params) {
  step("Analysing the NFP survey")
  out <- list(
    engagement            = engagement_summary(survey),
    engagement_dist       = engagement_distribution(survey),
    engagement_capacity   = engagement_by_capacity(survey),
    staff                 = staff_summary(survey),
    decisions_dist        = decisions_distribution(survey, params),
    decisions_headline    = decisions_headline(survey),
    value_choice          = value_choice_summary(survey),
    value_choice_capacity = value_choice_by_capacity(survey),
    integration           = integration_summary(survey, params),
    rankings              = rankings_summary(survey),
    free_text             = free_text_table(survey)
  )

  dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(out)) {
    write_csv(out[[nm]], file.path(params$output_dir, paste0("survey_", nm, ".csv")))
  }
  say(sprintf("wrote %d survey tables to output/", length(out)))
  out
}

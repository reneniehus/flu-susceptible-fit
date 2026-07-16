# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Load & tidy the de-identified NFP survey ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The export (data/survey_deidentified.xlsx, one sheet "Content") is a wide matrix
# with a two-row metadata banner, the question text on row 4, and one respondent
# per row below it. This loader turns that into two clean objects:
#   $codebook  - one row per analysed column: variable name, question, type, text
#   $responses - one tidy row per respondent, typed and with readable column names
# The survey went to EU/EEA National Focal Points (NFPs) for viral respiratory
# diseases; each row is one country's collective NFP response, de-identified.

# ---- |-the codebook: raw column index -> analysis variable ----
# `col` is the 1-based column position in the sheet; `type` drives the later casting.
#   engage  0-5 self-rated awareness/engagement with an ECDC modelling output
#   likert  Q7 decision-informing likelihood (Very unlikely ... Very likely)
#   agree   Q10 agreement (Strongly disagree ... Strongly agree)
#   staff   Q6 in-house modelling staff band
#   choice  Q8 which output is most valuable
#   rank    Q11/Q12 semicolon-ordered preference list
#   text    free-text
survey_codebook <- function() {
  tribble(
    ~col, ~variable,               ~question, ~type,     ~label,
    1L,   "q1_is_nfp",             "Q1",      "text",    "Confirmed NFP for viral respiratory diseases",
    5L,   "eng_covid_guidance",    "Q5.1",    "engage",  "ECDC COVID-19 risk assessments & technical guidance (2020-2023)",
    6L,   "eng_covid_forecast",    "Q5.2",    "engage",  "European COVID-19 Forecast Hub (2021-2024) - short-term forecasts",
    7L,   "eng_covid_scenario",    "Q5.3",    "engage",  "European COVID-19 Scenario Hub (2022-2023) - scenarios",
    8L,   "eng_respicast",         "Q5.4",    "engage",  "RespiCast (2023-current) - short-term forecasts",
    9L,   "eng_respicompass",      "Q5.5",    "engage",  "RespiCompass (2024-current) - scenario projections",
    10L,  "staff",                 "Q6",      "staff",   "In-house mathematical-modelling staff",
    12L,  "dec_healthcare",        "Q7",      "likert",  "Healthcare capacity planning",
    13L,  "dec_surveillance",      "Q7",      "likert",  "Inform surveillance activities",
    14L,  "dec_vaccination",       "Q7",      "likert",  "Planning vaccination campaigns",
    15L,  "dec_countermeasures",   "Q7",      "likert",  "Procurement of medical countermeasures (non-vaccine)",
    16L,  "dec_other",             "Q7",      "likert",  "Other activities",
    17L,  "dec_other_text",        "Q7b",     "text",    "Q7 'other', specified",
    18L,  "value_choice",          "Q8",      "choice",  "Most valuable output: RespiCast vs RespiCompass",
    19L,  "value_choice_text",     "Q8b",     "text",    "Reasons for the Q8 choice",
    20L,  "integration",           "Q10",     "agree",   "Clear mechanism to integrate modelling into decisions",
    21L,  "comms_forecast_rank",   "Q11",     "rank",    "Preferred channels for weekly short-term forecasts (RespiCast)",
    22L,  "comms_forecast_text",   "Q11b",    "text",    "Q11 'other', specified",
    23L,  "comms_scenario_rank",   "Q12",     "rank",    "Preferred channels for scenario insights (RespiCompass)",
    24L,  "comms_scenario_text",   "Q12b",    "text",    "Q12 'other', specified",
    25L,  "reflections",           "Q13",     "text",    "Additional reflections on ECDC respiratory-virus modelling"
  )
}

# ---- |-read one column out of the raw sheet, casting by its codebook type ----
# `raw` is the whole sheet read with col_names = FALSE; `data_rows` selects the
# respondent rows (everything below the question-text header row).
cast_survey_column <- function(raw, cb_row, data_rows, params) {
  x <- as.character(raw[[cb_row$col]][data_rows])
  x <- str_squish(x)                        # trim stray whitespace from the export
  x[x %in% c("", "-", "n/a", "N/A", "NA")] <- NA   # normalise the "blank" placeholders

  switch(cb_row$type,
    engage = suppressWarnings(as.integer(x)),                       # 0-5, higher = more engaged
    likert = factor(x, levels = params$levels_q7,    ordered = TRUE),
    agree  = factor(x, levels = params$levels_q10,   ordered = TRUE),
    staff  = factor(x, levels = params$levels_staff, ordered = TRUE),
    x                                                              # choice / rank / text stay character
  )
}

# ---- |-mother loader ----
load_survey <- function(params) {
  step("Loading the de-identified NFP survey")
  raw <- suppressMessages(read_excel(params$survey_xlsx, sheet = "Content",
                                     col_names = FALSE, .name_repair = "minimal"))

  cb        <- survey_codebook()
  data_rows <- (params$survey_header_row + 1L):nrow(raw)             # respondents sit below row 4

  # attach the verbatim question text from the header row (row 4) for the codebook
  cb$question_text <- vapply(cb$col, function(j) str_squish(as.character(raw[[j]][params$survey_header_row])),
                             character(1))

  # build the tidy respondent table, one cast column per codebook entry
  responses <- tibble(respondent_id = sprintf("r%02d", seq_along(data_rows)))
  for (i in seq_len(nrow(cb))) {
    responses[[cb$variable[i]]] <- cast_survey_column(raw, cb[i, ], data_rows, params)
  }

  # keep only genuine respondents (a confirmed NFP answer marks a real row)
  responses <- responses %>% filter(!is.na(q1_is_nfp))

  say(sprintf("respondents: %d (expected %d)", nrow(responses), params$n_expected_nfp))
  if (nrow(responses) != params$n_expected_nfp)
    warning("Survey respondent count differs from the expected de-identified export size.")

  list(codebook = cb, responses = responses)
}

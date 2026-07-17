# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Run configuration ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# One settings() function returning a params list, mirroring the "high-level user
# edits this one file" convention. Paths, the hub folder names, and the small set
# of naming rules that tell an ENSEMBLE / BASELINE apart from a genuine model all
# live here rather than being sprinkled through the analysis code.

settings <- function() {
  params <- list()

  # ---- |-Paths ----
  params$survey_xlsx <- here("data", "survey_deidentified.xlsx")   # the only committed input
  params$output_dir  <- here("output")                             # derived tables + artefact JSON
  params$figure_dir  <- here("output", "figures")                  # static figure companions

  # The two forecasting-hub clones live OUTSIDE the repo (they are ~1.7 GB of
  # submissions). Point here to wherever they were cloned. Default: a sibling
  # "hubs/" folder next to the repo; override with env var RESPICAST_HUBS_DIR.
  hubs_env <- Sys.getenv("RESPICAST_HUBS_DIR", unset = "")
  params$hubs_dir <- if (nzchar(hubs_env)) hubs_env else normalizePath(here("..", "hubs"), mustWork = FALSE)

  # ---- |-The forecasting hubs (current + archived) ----
  # name       : short id used in the tidy tables + artefact
  # folder     : the clone's directory name under hubs_dir
  # format     : "modern" = hubverse (model-output/<team>/<origin_date>-<team>.csv, target column)
  #              "legacy" = old EU COVID hub (data-processed/, compound "N wk ahead inc hosp" targets)
  # era        : "current" (live) or "archive" (superseded predecessor) -- for the handover story
  # indicators : the forecast target(s) that hub carries (verified against every file + git history)
  # Together these give a CONTINUOUS timeline: COVID hospitalisations since mid-2021 (legacy COVID hub
  # -> RespiCast-Covid19), and ILI/ARI since the 2023/24 season (flu/ari archives -> RespiCast-Syndromic).
  params$hubs <- tribble(
    ~name,            ~folder,                              ~format,  ~era,       ~indicators,
    "covid_archive",  "covid19-forecast-hub-europe_archive","legacy", "archive",  "COVID-19 cases / hospitalisations / deaths",
    "covid",          "RespiCast-Covid19",                  "modern", "current",  "COVID-19 hospitalisations",
    "flu_archive",    "flu-forecast-hub_archive",           "modern", "archive",  "ILI incidence",
    "ari_archive",    "ari-forecast-hub_archive",           "modern", "archive",  "ARI incidence",
    "syndromic",      "RespiCast-SyndromicIndicators",      "modern", "current",  "ILI + ARI incidence"
  )

  # ---- |-Model-role rules ----
  # A model-output/ (or data-processed/) folder is one team-model. A few are NOT ordinary
  # models and are counted apart from the "how many models" tally the survey cares about:
  #   ensemble : the hub's OFFICIAL combined product (team_model_designation == "ensemble").
  #              NB fjordhest-ensemble is a *participating team's* model (designation "primary"),
  #              not a hub ensemble, so it is deliberately NOT listed here -- it counts as a model.
  #   baseline : the reference baseline each hub ships to benchmark against.
  params$ensemble_models <- c("respicast-hubEnsemble", "EuroCOVIDhub-ensemble")
  params$baseline_models <- c("respicast-quantileBaseline", "EuroCOVIDhub-baseline")

  # ---- |-Human-readable target labels ----
  # modern hubs: the raw `target` string -> indicator label.
  # legacy COVID hub: the indicator is the last token of "N wk ahead inc|cum {case|hosp|death}",
  # mapped here from that token (see load_forecasts.R::read_legacy_index).
  params$target_labels <- c(
    "hospital admissions" = "COVID-19 hospitalisations",
    "ILI incidence"       = "ILI incidence",
    "ARI incidence"       = "ARI incidence",
    "case"                = "COVID-19 cases",
    "hosp"                = "COVID-19 hospitalisations",
    "death"               = "COVID-19 deaths"
  )

  # ---- |-Survey coding ----
  # The header row sits at row 4 of the export; respondents are the rows below it.
  params$survey_header_row <- 4L
  params$n_expected_nfp    <- 19L   # respondents in this de-identified export (sanity check)

  # Ordered factor levels for the Likert-style questions, so plots and summaries
  # sort them meaningfully rather than alphabetically.
  params$levels_q7    <- c("Very unlikely", "Unlikely", "Unsure", "Likely", "Very likely")
  params$levels_q10   <- c("Strongly disagree", "Disagree", "Neither agree nor disagree",
                           "Agree", "Strongly agree")
  params$levels_staff <- c("0 staff", "1-5 staff", "6-10 staff", ">10 staff")

  return(params)
}

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

  # ---- |-The two hubs ----
  # name        : short id used in the tidy tables + artefact
  # folder      : the clone's directory name under hubs_dir
  # indicators  : the forecast target(s) that hub carries (verified against every file + git history)
  params$hubs <- tribble(
    ~name,        ~folder,                         ~indicators,
    "covid",      "RespiCast-Covid19",             "COVID-19 hospitalisations",
    "syndromic",  "RespiCast-SyndromicIndicators", "ILI + ARI incidence"
  )

  # ---- |-Model-role rules ----
  # A hub folder under model-output/ is one team-model. Three of those folders are
  # NOT ordinary models and must be counted separately from the "how many models"
  # tally the survey cares about:
  #   ensemble : combined multi-model product (the hub's flagship external output)
  #   baseline : the reference quantile baseline every hub ships to benchmark against
  params$ensemble_models <- c("respicast-hubEnsemble", "fjordhest-ensemble")
  params$baseline_models <- c("respicast-quantileBaseline")

  # ---- |-Human-readable target labels ----
  # the raw `target` strings in the files -> the indicator labels used everywhere else
  params$target_labels <- c(
    "hospital admissions" = "COVID-19 hospitalisations",
    "ILI incidence"       = "ILI incidence",
    "ARI incidence"       = "ARI incidence"
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

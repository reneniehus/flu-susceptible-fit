# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 00_main.R -- build every analysis table + the artefact data, end to end ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run from the repo root:  Rscript code/00_main.R
#
# Two independent analyses, one shared session:
#   A. the de-identified NFP survey  (code/02_survey/)  -> output/survey_*.csv
#   B. the RespiCast forecast hubs   (code/03_hubs/)    -> output/hub_*.csv
# Both are then folded into one artefact JSON (code/04_export/) that the HTML
# dashboards embed. See README.md for the full picture and how to get the hub data.

gc()

# ---- |-Set up ----
source("code/01_support/setup.R")
source("code/01_support/config.R"); params <- settings()

# ---- |-Source the analysis modules ----
source("code/02_survey/load_survey.R")      # read + tidy the survey
source("code/02_survey/analyse_survey.R")   # survey summary tables
source("code/03_hubs/load_forecasts.R")     # scan the hub submissions
source("code/03_hubs/analyse_coverage.R")   # weeks x indicator x #models coverage
source("code/04_export/build_artefact_data.R")

# ---- |-A. Survey: modelling capacity & the value of ECDC forecasting ----
survey     <- load_survey(params)
survey_out <- run_survey_analysis(survey, params)

# ---- |-B. Hubs: two seasons of forecast coverage ----
submissions  <- load_hub_forecasts(params)
coverage_out <- run_coverage_analysis(submissions, params)

# ---- |-Assemble the artefact data (one JSON the dashboard embeds) ----
# The qualitative theme synthesis is produced by the verification/coding workflow and
# cached to output/themes.json; fold it in when present, otherwise build without it.
themes_path <- file.path(params$output_dir, "themes.json")
themes <- if (file.exists(themes_path)) jsonlite::read_json(themes_path, simplifyVector = FALSE) else NULL
build_artefact_data(survey, survey_out, coverage_out, params, themes = themes)

step("Done")
say("survey tables + hub tables -> output/*.csv")
say("artefact data            -> output/artefact_data.json")
say("build the figures with:   Rscript code/05_figures/fig_survey.R  &  fig_coverage.R")

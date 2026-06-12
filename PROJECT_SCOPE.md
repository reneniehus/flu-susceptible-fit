---
editor_options: 
  markdown: 
    wrap: 72
---

# Project Scope

A clone-and-go R template for loading ERVISS / RespiCompass
respiratory-virus surveillance data (plus contact, vaccine and
demography data) and turning it into tidy, model-ready tables with a
data-quality / dynamics report. The data layer is the durable core;
modelling is a deliberately light scaffold to graduate into production
per project.

See `documentation/quickstart.md` (how to run) and
`documentation/data_overview.md` (what data is present).

## Active production path

1.  `code/00_main.R` — entrypoint / orchestrator
2.  `code/02_settings/settings_version0.R` — runtime settings (`params`)
3.  `code/01_main_supporting/setup.R` — libraries + shared helpers
4.  `code/01_main_supporting/validate.R` — input-data contracts
5.  `code/01_main_supporting/load_data.R` — data loading and caching
    (`data`)
6.  `code/01_main_supporting/gen_model_input.R` — canonical
    long/wide/season-summary tables (`models_in`)
7.  `code/01_main_supporting/eyeballing.R` +
    `code/03_report/eyeballing_report.Rmd` — quality/dynamics report

## In scope

-   Loading and caching multi-pathogen ERVISS/RespiCompass streams
    (registry-driven, extendable).
-   Tidy canonical tables: `data_timeseries_long`,
    `data_timeseries_wide`, `data_season_summary`.
-   Multi-pathogen ILI+ (Influenza / SARS-CoV-2 / RSV) and data-quality
    measures.
-   Contact, vaccine and demography data (from committed snapshots).

## Out of scope (for production now)

-   Heavy modelling. `run_model.R` is an empty orchestration stub and
    `model_*.R` are parked single-model templates (SIR, ARIMA,
    last-year-burden) to be wired in when needed.
-   Hub submission / reporting beyond the eyeballing report
    (`process_and_save.R`, `send_report.R` are scaffolds).

## Reproducibility & tests

-   `renv.lock` pins all dependencies; `renv::restore()` reproduces the
    environment.
-   `tests/testthat/` checks the data contracts and key invariants; run
    `Rscript run_tests.R`.

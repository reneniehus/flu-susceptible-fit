# Data overview

Two datasets, one committed and one external.

## 1. The NFP survey — `data/survey_deidentified.xlsx`

One sheet ("Content"). The export has a two-row metadata banner, the **question text on row 4**, and
**one respondent per row below it** (rows 5 down). This project reads **19 respondents**. Each row is
one country's collective NFP response, de-identified — there is no per-respondent country label.

`code/02_survey/load_survey.R` turns the wide matrix into a tidy respondent table using an explicit
**codebook** (raw column index → analysis variable). The analysed columns:

| Q | Variable | Type | Meaning |
|---|---|---|---|
| Q1 | `q1_is_nfp` | text | confirmation of NFP role (all "Yes") |
| Q5.1 | `eng_covid_guidance` | 0–5 | engagement: ECDC COVID-19 risk assessments & guidance (2020–2023) |
| Q5.2 | `eng_covid_forecast` | 0–5 | engagement: European COVID-19 Forecast Hub (2021–2024) — *forecasting* |
| Q5.3 | `eng_covid_scenario` | 0–5 | engagement: European COVID-19 Scenario Hub (2022–2023) — *scenario* |
| Q5.4 | `eng_respicast` | 0–5 | engagement: RespiCast (2023–) — *forecasting* |
| Q5.5 | `eng_respicompass` | 0–5 | engagement: RespiCompass (2024–) — *scenario* |
| Q6 | `staff` | band | in-house modelling staff: `0 staff` / `1-5 staff` / `>10 staff` |
| Q7 | `dec_*` | Likert | likelihood modelling informs 5 actions (Very unlikely … Very likely) |
| Q8 | `value_choice` | choice | most valuable: RespiCast / RespiCompass / Both / Neither |
| Q10 | `integration` | Likert | agreement there is a clear mechanism to integrate modelling |
| Q11/Q12 | `comms_*_rank` | ranking | preferred channels for forecasts / scenarios (";"-ordered) |
| Q8b, Q13, … | `*_text` | free text | open responses |

**The Q5 0–5 scale** is treated as an **ordinal awareness/engagement** score (0 = not aware / no use,
higher = more engaged). The exact anchor wording is not present in the de-identified export, so only
the ordering and distribution are interpreted, never a conversion to "% who find it useful".

## 2. The forecasting hubs — external clones under `../hubs/`

**Five** hubs, chosen so the record is *continuous* — each live RespiCast hub is paired with the
archived predecessor it replaced. Two file formats (dispatched on `format` in `config.R$hubs`):

| Hub | Format | Era | Indicator(s) | Window (ISO weeks) |
|---|---|---|---|---|
| covid19-forecast-hub-europe_archive | legacy | archive | COVID-19 cases / hospitalisations / deaths | 2021-02 → 2024-10 |
| RespiCast-Covid19 | modern | current | COVID-19 hospitalisations | 2024-10 → 2026-06 |
| flu-forecast-hub_archive | modern | archive | ILI incidence | 2023-11 → 2024-05 |
| ari-forecast-hub_archive | modern | archive | ARI incidence | 2023-12 → 2024-05 |
| RespiCast-SyndromicIndicators | modern | current | ILI + ARI incidence | 2024-10 → 2026-07 |

**Modern (hubverse) format** — `model-output/<team>/<origin_date>-<team>.csv`, with an `origin_date`
(Wednesday round) and a `target` column naming the indicator. **Legacy (old EU COVID hub) format** —
`data-processed/<team>/<forecast_date>-<team>.csv`, with a `forecast_date` (Monday round) and compound
targets like `"2 wk ahead inc hosp"`; the indicator is the last token (`case` / `hosp` / `death`).

Three quirks, all handled in `code/03_hubs/load_forecasts.R` and found by inspecting the raw files:
1. **Column order is not constant across modern files** — columns are parsed **by name**, never by position.
2. A handful of submissions are **header-only** (no data rows) — skipped and counted.
3. Legacy rounds fall on **Mondays**, modern rounds on **Wednesdays** — both are snapped to the
   **ISO-week Monday** (`week_monday()`) so the two eras land on one weekly grid; that shared grid is
   what makes gap-detection (the "uninterrupted since 2021" test) meaningful.

There is **no "COVID cases" target** in the RespiCast era — it lived only in the archived hub. The
syndromic hub's git history shows it has **only ever** carried ILI and ARI.

### Model roles

A `model-output/` (or `data-processed/`) folder is one team-model. A few are **not** ordinary models
and are counted apart:

- **ensembles (the official hub product):** `respicast-hubEnsemble` (RespiCast), `EuroCOVIDhub-ensemble`
  (legacy). **Note:** `fjordhest-ensemble` is a *participating team's* model (Fjordhest, NIPH;
  `team_model_designation: primary`), not a hub ensemble — so it counts as a **model**, not an ensemble.
- **baseline:** `respicast-quantileBaseline`, `EuroCOVIDhub-baseline`.

Everything the analysis calls "models" excludes the ensembles and baselines.

### The other org repositories (assessed, and why they are / aren't used)

The `european-modelling-hubs` org has 23 repositories. Beyond the five above, the rest are **not**
sources of time-stamped forecast submissions and are out of scope for a coverage analysis: the COVID/
RespiCompass **scenario** hubs (scenario projections, not weekly forecasts — `RespiCompass`,
`covid19-scenario-hub-europe-1`), **websites/viz** (`*-website`, `*-viz`, `actions-dashboard`),
**tooling / validation / R packages** (`hub-tools`, `HubValidations`, `EuroForecastHub`,
`*-baseline`, `*-validations`, `HubSubmissionApp`, `modelling_setup`, `respicompass-resources`,
`.github`), and **per-team auto-submission** repos (`*-submissions`, `autosubmission-*`). The scenario
hubs are related to the *survey* (Q5.3 / Q5.5) but hold scenarios, not the weekly forecasts this
delivery analysis counts.

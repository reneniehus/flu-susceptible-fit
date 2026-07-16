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

## 2. The RespiCast hubs — external clones under `../hubs/`

Both follow the [hubverse](https://hubverse.io) layout:
`model-output/<team-model>/<origin_date>-<team-model>.csv`, one file per team per weekly round.

| Column | Meaning |
|---|---|
| `origin_date` | the weekly forecast round (a Wednesday); the analysis unit |
| `target` | the indicator (see below) |
| `target_end_date`, `horizon` | 1–4 weeks ahead |
| `location` | EU/EEA country (ISO2) |
| `output_type`, `output_type_id`, `value` | median + quantile forecast values (not used here) |

Two quirks, both handled in `code/03_hubs/load_forecasts.R` and found by inspecting the raw files:
1. **Column order is not constant across files** — columns are always parsed **by name**, never by
   position.
2. A handful of submissions are **header-only** (no data rows) — skipped and counted.

### Targets (verified against every file and the git history)

| Hub | Target(s) | Window |
|---|---|---|
| RespiCast-Covid19 | `hospital admissions` → **COVID-19 hospitalisations** | 2024-10-23 → 2026-06-24 (88 rounds) |
| RespiCast-SyndromicIndicators | `ILI incidence`, `ARI incidence` | 2024-10-23 → 2026-07-15 (91 rounds) |

There is **no "COVID cases" target** in either hub. COVID-19 hospitalisations sit in their own hub;
the syndromic hub's git history shows it has **only ever** carried ILI and ARI. (The COVID hub's
config history still contains an "ILI incidence" template stanza — a trace of the RespiCast
reorganisation that split the indicators into separate repositories.)

### Model roles

A `model-output/` folder is one team-model. Three are **not** ordinary models and are counted apart:

- **ensembles:** `respicast-hubEnsemble` (both hubs), `fjordhest-ensemble` (syndromic) — the combined
  products external stakeholders actually see.
- **baseline:** `respicast-quantileBaseline` — the reference every hub ships to benchmark against.

Everything the analysis calls "models" excludes these three.

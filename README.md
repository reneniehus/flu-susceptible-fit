# respicast-value — demand & delivery of ECDC respiratory-virus forecasting

Two views of the same system, built from source and rendered into one interactive dashboard:

- **Demand** — a survey of **19 EU/EEA National Focal Points (NFPs)** for viral respiratory diseases,
  on their in-house modelling capacity and the value they place on ECDC's short-term forecasts
  (**RespiCast**) and seasonal scenarios (**RespiCompass**).
- **Delivery** — **five** European forecasting hubs reconstructed submission-by-submission, back to 2021:
  the archived **EU COVID-19 Forecast Hub** and its **RespiCast** successors, plus the archived **flu**
  and **ARI** hubs — which indicators were forecast, in which weeks, by how many models and the ensemble.

The brief was to focus on the **value (and potential value) of forecasting / nowcasting by ECDC for
external stakeholders** such as national public health institutes — so the whole analysis circles one
axis: national modelling **capacity**, because that is what decides whether an ECDC forecast is a
convenience or the only forecast a country has.

**Interactive dashboard:** https://claude.ai/code/artifact/83633670-7442-4950-b097-d6556709f5e9

## Headline findings

**Demand (survey).**
- **13 of 19** focal points have **no in-house mathematical modeller at all**; only one institute has a
  team of any size. For most countries "using a forecast" can only mean using someone else's.
- Engagement with the newer flagship products is modest — RespiCast mean **2.11/5**, RespiCompass
  **1.35/5** — but they have had far less time in market than the older COVID-era outputs (a recency
  confound, not a value verdict).
- When forced to choose, focal points prefer **RespiCast (forecasts) over RespiCompass (scenarios)
  5:1**; a majority expect **both** to be useful (an *expectation* — many have not yet used them).
- Engagement is **flat across capacity bands** — the institutes that would gain most engage no more
  than better-resourced ones. The value is **potential, not yet realised**; the binding constraint is
  demand-side (awareness, interpretation, a route into the decision), not supply.
- Forecasting is expected to inform **surveillance** most (58% likely) and **healthcare-capacity
  planning** — the canonical forecasting use — least (16%). Only **~half** report a clear mechanism to
  integrate modelling into decisions.

**Delivery (five hubs, 2021–2026).**
- **The claim tested** — *"COVID-19 hospitalisation forecasts have been produced almost uninterrupted
  since 2021, first by the COVID-19 hub and then within RespiCast"*: **supported, with one refinement.**
  Hospitalisation forecasts ran from **26 Jul 2021 → 22 Jun 2026**, covering **255 of 257 ISO weeks
  (99.2%)** with only **two isolated one-week gaps** and a **seamless** archive→RespiCast handover
  (consecutive weeks). The refinement: "since 2021" is precisely **late July 2021** — COVID *cases and
  deaths* were forecast from February 2021, but hospitalisations began that summer. And the record is
  deep only recently: a **median of ~4 models** per week in the EU-COVID-hub era, roughly doubling to
  ~8.5 under RespiCast.
- **ECDC narrowed COVID forecasting.** The archived EU hub forecast **cases, hospitalisations and
  deaths** — cases/deaths with the deepest fields (up to **32 models**, 32 countries) — but those two
  **ended at the October-2024 handover**; RespiCast carried only hospitalisations forward.
- **Syndromic indicators** (ILI, ARI) reach back one further season: the archived flu and ARI hubs ran
  the **2023/24** winter before RespiCast-SyndromicIndicators. Both run winter-only, with a summer-long
  gap between eras. ILI still draws ~2× the models of ARI from the same files.
- No **"COVID cases"** target exists in the RespiCast era — that indicator lived only in the archived hub.

> Every number above is regenerated from source by `code/00_main.R`. Key figures were independently
> re-derived from the raw files by a verification pass (see `documentation/reflections.md`); the survey
> counts, target lists, date/round counts, country counts, per-week model counts **and the
> hospitalisation continuity** all matched.

## Quick start

```r
# 1. get the five hub repositories (NOT vendored here -- ~4 GB of submissions).
#    clone them next to this repo, into ../hubs/  (or set RESPICAST_HUBS_DIR). Shallow (--depth 1)
#    is fine: the current tree already holds every historical submission.
#      git clone --depth 1 https://github.com/european-modelling-hubs/covid19-forecast-hub-europe_archive.git ../hubs/covid19-forecast-hub-europe_archive
#      git clone --depth 1 https://github.com/european-modelling-hubs/RespiCast-Covid19.git                    ../hubs/RespiCast-Covid19
#      git clone --depth 1 https://github.com/european-modelling-hubs/RespiCast-SyndromicIndicators.git        ../hubs/RespiCast-SyndromicIndicators
#      git clone --depth 1 https://github.com/european-modelling-hubs/flu-forecast-hub_archive.git             ../hubs/flu-forecast-hub_archive
#      git clone --depth 1 https://github.com/european-modelling-hubs/ari-forecast-hub_archive.git             ../hubs/ari-forecast-hub_archive

# 2. build every analysis table + the artefact data
Rscript code/00_main.R

# 3. (re)build the dashboard page and the static figures
Rscript code/05_artefact/build_pages.R      # -> artefact/dashboard.html
Rscript code/05_figures/fig_survey.R        # -> output/figures/survey.png
Rscript code/05_figures/fig_coverage.R      # -> output/figures/coverage.png
```

The committed `output/*.csv` + `output/artefact_data.json` already hold the derived results, so the
dashboard renders without re-cloning the hubs; step 1 is only needed to regenerate them from scratch.

### Reproducing the interactive dashboard

The published dashboard (`artefact/dashboard.html`) is fully reproducible from committed source:

```sh
./reproduce.sh          # full: data -> artefact JSON -> dashboard + figures (needs the 5 hub clones)
./reproduce.sh page     # just rebuild artefact/dashboard.html from the committed JSON (no hubs needed)
```

`code/05_artefact/build_pages.R` injects `output/artefact_data.json` into
`code/05_artefact/dashboard_template.html` to produce `artefact/dashboard.html` — a single,
self-contained file (all data embedded, all charts hand-built inline, no external requests). `./reproduce.sh page`
regenerates it **byte-identically** from the committed JSON, so the artefact can always be rebuilt and
re-published without the hub data. Publish it as-is (e.g. open in a browser, or re-publish as an Artifact).

## Layout

```
data/survey_deidentified.xlsx   the ONE committed input (de-identified NFP survey export)
code/00_main.R                  orchestrator: survey + hubs -> output/ tables + artefact JSON
code/01_support/                setup.R (libraries, palette), config.R (settings/params)
code/02_survey/                 load_survey.R (tidy + codebook), analyse_survey.R (summary tables)
code/03_hubs/                   load_forecasts.R (scan submissions), analyse_coverage.R (coverage)
code/04_export/                 build_artefact_data.R (tables -> one JSON blob)
code/05_artefact/               dashboard_template.html (+ build_pages.R injects the data)
code/05_figures/                fig_survey.R, fig_coverage.R (PNG companions)
output/                         derived tables (survey_*.csv, hub_*.csv), artefact_data.json, figures/
artefact/dashboard.html         the self-contained interactive dashboard (published as an Artifact)
documentation/                  data_overview, decisions, findings, reflections (see table below)
```

## Documentation

| File | Holds |
|---|---|
| `README.md` | what this is, headline findings, how to run, layout |
| `PROJECT_SCOPE.md` | aim, in/out of scope, data sources |
| `documentation/data_overview.md` | the two datasets: survey structure/codebook and hub schema |
| `documentation/findings.md` | the results in full (survey + coverage), with caveats |
| `documentation/decisions.md` | **why** — key analysis / coding / design decisions (append-only) |
| `documentation/reflections.md` | interpretation, the verification pass, overclaims to avoid, next analyses |
| inline header comments | why each script / function works the way it does |

## Reproducibility & caveats

Public data throughout. The survey has **19 respondents** (of ~30 EU/EEA NFPs) — every share moves ~5
points per person, so read counts, not precise percentages. Q5's 0–5 scale is read **ordinally**
(higher = more engaged); the exact anchor wording is not in the de-identified export. Hub coverage
counts a submission's **presence**, not its skill or accuracy. See `documentation/reflections.md` for
the full list of caveats and the results of the independent verification pass.

# Using real data instead of the simulation

The analysis tool box never touches the simulation object. It consumes a small set of **tidy data
frames in a fixed schema**, produced for a simulated run by `as_analysis_input(sim)`. To analyse
**real** surveillance data, build a list with the same frames and column names and hand it to the same
tools — they cannot tell the difference. What you lose with real data is only the `score_*()` layer,
because the truth is unknown (which is the whole reason the playground exists).

## The schema

`as_analysis_input()` returns a named list. Each analysis tool uses only the frames it needs, so you
only have to supply the frames for the tools you want to run.

| Element | Columns | Used by |
|---|---|---|
| `cases_by_onset` | `location, day, date, cases` (right-truncated by onset date) | growth, Rt, forecast |
| `cases_by_report` | `location, day, date, cases` (by report date; complete) | growth, CFR |
| `reporting_triangle` | `location, onset_day, report_day, cases` | nowcast, forecast |
| `deaths` | `location, day, date, deaths_by_date, deaths_by_onset` | CFR |
| `admissions` | `location, day, date, admissions` | forecast |
| `detected_imports` | `country, day, date, detected_imports, surveillance_quality` | catchment, importation risk |
| `variant_cases` | `location, day, date, sequenced, variant` | variant selection |
| `flight_volumes` | `country, day, date, volume` (observed) | catchment, importation risk |
| `delays` | a list of `epidist` objects: `generation_interval, incubation, onset_to_death, onset_to_report, onset_to_admission` | most tools |
| `as_of` | integer day index of the data cutoff | nowcast, CFR, forecast |
| `countries`, `source_code` | character | bookkeeping |

Conventions: `day` is an integer index (0-based is fine; tools only use ordering and differences);
`location`/`country` are free-text codes; counts are non-negative integers. The delay distributions are
the analyst's **assumed** (literature) delays — build them with `epidist_gamma()` /
`epidist_lognormal()` or, if `{epiparameter}` is installed, `as_epidist(<epiparameter object>)`.

## Minimal example (real data, no simulator)

```r
source("setup.R")

# your own surveillance tables, in the schema
real <- list(
  cases_by_report = data.frame(location = "MYCOUNTRY", day = 0:60, date = as.Date("2027-03-01") + 0:60,
                               cases = my_reported_cases),
  deaths          = data.frame(location = "MYCOUNTRY", day = 0:60, date = as.Date("2027-03-01") + 0:60,
                               deaths_by_date = my_deaths, deaths_by_onset = NA),
  delays          = list(generation_interval = epidist_gamma("gi", 5.2, 1.7),
                         onset_to_death      = epidist_gamma("o2d", 15, 6.6)),
  as_of           = 60, source_code = "MYCOUNTRY"
)

growth_analysis(real, "MYCOUNTRY", window = 0:25)     # Phase 0 growth -> R
cfr_static(real, "MYCOUNTRY", as_of = 60)             # Phase 1 delay-adjusted CFR
```

`test-real-data.R` checks exactly this: a tool gives the same answer whether fed the sim adapter or a
hand-built schema list, and it recovers a planted growth rate / selection coefficient from purely
hand-made data.

## Substitution points, per data source

- **Case line lists / aggregates → `cases_by_report`, `cases_by_onset`, `reporting_triangle`.** If you
  have a line list with onset and report dates, cross-tabulate onset × report into the triangle; the
  by-onset and by-report frames are its margins. If you only have report-date counts, supply
  `cases_by_report` and skip the nowcast.
- **Death registrations → `deaths`.** `deaths_by_date` is the daily death count; `deaths_by_onset` (the
  onset day of fatal cases) is optional and only used by the onset-based CFR path.
- **Hospital admissions → `admissions`.** Daily admission counts; needed only for the capacity forecast.
- **Traveller screening / imported-case counts → `detected_imports`.** Include a
  `surveillance_quality` covariate (0–1) per country if you have one — the importation-risk regression
  anchors on the high-quality countries.
- **Flight statistics → `flight_volumes`.** Daily (or interpolated-to-daily) outbound volume from the
  source to each country.
- **Genomic surveillance → `variant_cases`.** `sequenced` = samples typed that day, `variant` = of those,
  how many were the variant.
- **Delays → `delays`.** From the literature or `{epiparameter}`; pass a different `epidist` to any tool
  to test sensitivity to delay misspecification.

# Design decisions & rationale

The README records *what* the project does; this file records *why* the key analysis, coding and design
choices were made, so the reasoning survives beyond commit messages. Append new decisions as they are
made — each entry: the **decision**, the reason, and the main alternative considered.

## Framing

- **Organise everything around one axis: in-house modelling capacity (Q6).** The brief is the value of
  ECDC forecasting to external stakeholders; capacity is the variable that decides whether a forecast
  is a convenience or the only forecasting a country has. *Alternative:* report each question
  independently — rejected as a pile of facts with no through-line.

- **Read the two datasets as demand and supply of one system.** The survey's RespiCast/RespiCompass
  questions ask NFPs to value the exact products the hubs deliver, so pairing them lets "is it valued?"
  meet "is it delivered?". *Alternative:* two unrelated analyses — rejected; the bridge is the insight.

## Survey coding

- **CORRECTION — Q5 measures AWARENESS, not engagement/use.** After the official question wording was
  supplied — *"How aware are you of respiratory virus burden modelling work done at ECDC?"* (0–5 slider,
  0 = not aware, 5 = fully aware) — every "engagement / actively used / uptake" reading of the Q5 scores
  was corrected to **awareness** across all deliverables (dashboard, decision note, supplement, figures,
  docs). The internal variable/key names (`eng_*`, `engagement*`) are kept to avoid destabilising the
  reproducible pipeline, but are documented as historical at their definitions. Awareness is a *necessary,
  not sufficient*, condition for value; the survey does **not** measure actual use (only the open text does).

- **CORRECTION — the instrument is a general modelling survey, not a RespiCompass evaluation.** The
  official intro shows it was run by the ECDC Respiratory Viruses & Legionella group + Modelling team to
  learn whether/how countries use modelling and **which ECDC outputs they would prefer** (results to a
  network meeting, October 2024). So Q8 (forecasts vs scenarios preference) is an **in-scope** question,
  not a by-product — the earlier "drawn from a RespiCompass-designed instrument, not relied on" caveat was
  wrong and has been reversed.

- **Read the Q5 awareness gradient as recency.** Awareness falls from guidance (since 2020) to RespiCast
  (2023) to RespiCompass (2024) because those are launch dates; lower awareness of newer products is
  mostly *time in circulation*, not lesser value or use. *Alternative:* present it as a value ranking —
  rejected as misleading.

- **The awareness gap is concentrated where capacity is absent; do not extend it to "use".** The Q5×Q6
  cross-tab shows awareness roughly flat across staff bands (2.4 / 2.3 / 3.0, the last n=1), and all three
  institutes entirely unaware of RespiCast have no team — an awareness/distribution gap, cheap to close.
  It says nothing about whether awareness converts into use. *Alternative (an earlier draft):* reading the
  scores as "engagement/use" and inferring demand conversion — **retracted** as a misread of Q5.

## Hub coverage

- **Parse hub CSVs by column name, not position.** Column order is not constant across files, so
  positional reads silently mixed up `target`/`location`/`horizon`. *Alternative:* trust the header
  order — rejected after finding the inconsistency.

- **Count ensembles and the baseline separately from "models".** The survey cares about "how many
  models (+ ensemble)"; the ensemble and the reference quantile-baseline are products, not contributing
  models, and are tallied apart (`role` in the loader). *Alternative:* count every folder equally —
  rejected as it would inflate the model count by the hub's own outputs.

- **Coverage = presence, and say so everywhere.** We measure whether a submission exists, not its skill.
  Every coverage view carries this caveat so "covered" is never misread as "accurate". *Alternative:*
  bring in scores now — deferred (needs the quantiles + truth data; the highest-value next step).

- **Include the archived predecessor hubs; pair each live hub with the one it replaced.** To test the
  "COVID-19 hospitalisation forecasts since 2021" claim we added the archived EU COVID hub
  (`covid19-forecast-hub-europe_archive`), and — for symmetry and to extend ILI/ARI back a season — the
  archived flu and ARI hubs. The other 18 org repos are scenario hubs, websites, tooling or
  auto-submission repos with no time-stamped forecast submissions, so they are out of scope for a
  coverage analysis (documented in `data_overview.md`). *Alternative:* the COVID archive alone —
  rejected; the flu/ARI archives were cheap (modern format, ~120 MB) and complete the syndromic timeline.

- **Snap every round to its ISO-week Monday, then detect gaps on that shared grid.** Legacy COVID-hub
  rounds fall on Mondays, modern RespiCast rounds on Wednesdays; without a common key the archive→live
  handover looks like a 9-day gap when it is actually two consecutive ISO weeks. `week_monday()` maps
  both to the week's Monday; continuity is then "is every Monday from first to last present?".
  *Alternative:* use raw dates — rejected; it fabricates handover gaps and can't align the two formats.

- **Count `fjordhest-ensemble` as a model, not an ensemble.** Its metadata is `team_name: Fjordhest`
  (Norwegian Institute of Public Health), `team_model_designation: primary` — a participating team that
  happens to use an ensemble method, not the hub's official product. Only the hub-designated ensembles
  (`respicast-hubEnsemble`, `EuroCOVIDhub-ensemble`) count as "the ensemble". This corrects the first
  pass, which had wrongly listed fjordhest as an ensemble (it inflated "ensemble present" and undercounted
  models by one in the weeks it submitted).

- **Parse the legacy compound target by its last token.** Old EU-hub targets encode indicator + horizon
  in one string (`"2 wk ahead inc hosp"`). We extract the indicator from the trailing `case|hosp|death`
  token and the horizon from the leading integer, then map to the same indicator labels the modern hubs
  use — so both formats flow into one `submissions` table. *Alternative:* a separate legacy schema —
  rejected; one shared schema keeps every downstream view format-agnostic.

## Artefact & design

- **One unified dashboard, two chapters, sticky nav — not two separate pages.** The strongest telling
  is demand → bridge → delivery in a single shareable link with one design system. *Alternative:* two
  artefacts — rejected as chrome duplication that hides the bridge (the point).

- **A "forecasting-desk" identity: humanist sans for prose, monospace for all data / week-codes /
  counts.** Grounds the look in the subject (a surveillance console) rather than a generic template, and
  the mono/sans split is a genuine, available-font pairing. Indicator hues (COVID blue, ILI orange, ARI
  violet) are reserved for data and never for chrome.

- **Charts hand-built as inline SVG, no libraries; validated palette; both themes.** The Artifact
  sandbox blocks external scripts/fonts, so everything is self-contained and the data is embedded as one
  JSON blob. The categorical (indicator) palette passes the design-system validator in light and dark;
  the page re-renders on theme change so computed colours stay correct. *Alternative:* a charting CDN —
  impossible under the CSP.

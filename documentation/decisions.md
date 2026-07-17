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

- **Treat Q5's 0–5 as ordinal engagement, never as a usefulness percentage.** The de-identified export
  does not carry the scale's anchor wording, so we interpret ordering and distribution only. This also
  guards against the recency confound (below). *Alternative:* map to "% who find it useful" — rejected
  as unsupported by the data.

- **Flag the recency confound on the engagement gradient explicitly.** Engagement falls from guidance
  (since 2020) to RespiCast (2023) to RespiCompass (2024); because those are also their launch dates,
  lower engagement with newer products is largely *time in market*, not lesser value. Every place the
  gradient appears carries this caveat. *Alternative:* present it as a value ranking — rejected as
  misleading.

- **Do not claim capacity-poor NFPs engage more.** The Q5×Q6 cross-tab shows short-term-forecasting
  engagement roughly flat across staff bands (2.4 / 2.3 / 3.0, the last n=1). The honest value story is
  the *need* gap plus stated preference, not differential engagement — so the framing is "near-universal
  need, flat uptake → demand-side bottleneck". This correction came directly from the verification/audit
  pass. *Alternative (an earlier draft):* "engagement rises where capacity is thin" — **retracted** as
  unsupported.

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

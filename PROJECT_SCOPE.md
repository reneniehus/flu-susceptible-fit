# Project scope — the value of ECDC forecasting to national focal points

## Aim

Characterise the **value, and potential value, of ECDC forecasting and nowcasting exercises for
external stakeholders** — principally EU/EEA national public health institutes — as ECDC develops
**RespiCast** (short-term forecasts) and **RespiCompass** (scenario modelling). The analysis pairs two
independent evidence bases:

1. **Demand.** A de-identified survey of National Focal Points (NFPs) for viral respiratory diseases,
   probing national modelling capacity, awareness of and engagement with ECDC modelling outputs, the
   decisions those outputs might inform, and the relative value NFPs place on forecasts vs scenarios.
2. **Delivery.** The public RespiCast forecasting-hub repositories, which record every model
   submission and every ensemble — used here to measure coverage: which indicators were forecast, in
   which weeks, and by how many models.

## Research questions

- **RQ1 (capacity).** How much modelling can EU/EEA NFPs do in-house, and how does that shape their
  reliance on ECDC outputs?
- **RQ2 (value).** How much do NFPs engage with, and expect to value, ECDC forecasting (RespiCast) vs
  scenarios (RespiCompass)? Which national decisions might they inform, and what stops them being used?
- **RQ3 (delivery).** Over the last two seasons, what did the RespiCast hubs actually produce — by
  indicator, week, model count and ensemble — and how well does that supply match the demand in RQ1–2?

## In scope

- Descriptive analysis of the de-identified NFP survey (closed and open responses).
- Reconstruction of forecast **coverage** across the two RespiCast hubs (COVID-19 hospitalisations;
  ILI and ARI incidence): weeks × indicators × number of models (+ ensemble, + baseline), model
  participation over time, and country reach.
- A single interactive dashboard artefact bridging demand and delivery, plus static figure companions.
- Independent verification of the headline numbers against the raw source files.

## Out of scope

- **Forecast skill / accuracy / calibration.** Coverage counts the *presence* of a submission, not how
  good it was. Nothing here evaluates whether the ensemble beat the baseline or was well-calibrated
  (this is the single most valuable follow-up — see `documentation/reflections.md`).
- Re-identification of survey respondents or countries (the export is de-identified by design).
- Causal inference. With 19 respondents (~30 EU/EEA NFPs) the survey is descriptive and
  hypothesis-generating; no significance testing or confidence intervals are reported.

## Data sources

- **Survey:** `data/survey_deidentified.xlsx` — ECDC RespiCompass 2024 NFP survey, de-identified. The
  only committed input.
- **Hubs (external, cloned separately):**
  - [RespiCast-Covid19](https://github.com/european-modelling-hubs/RespiCast-Covid19) — COVID-19
    hospital admissions.
  - [RespiCast-SyndromicIndicators](https://github.com/european-modelling-hubs/RespiCast-SyndromicIndicators)
    — ILI and ARI incidence.

Both hubs are part of the RespiCast consortium (ECDC with the ISI Foundation and LSHTM). Their forecast
targets, the RespiCast (short-term) and RespiCompass (scenario) programmes, are the exact outputs the
survey asks NFPs to value — which is what lets the two datasets be read as one story.

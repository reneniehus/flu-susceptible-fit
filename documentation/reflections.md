# Reflections — verification, caveats, and where the value really is

## The strongest, defensible value story

**ECDC forecasting is best understood as capacity substitution — a supranational public good that most
member states cannot produce for themselves — whose value is currently capped by uptake, not supply.**

The two datasets interlock. *Demand:* 68% of responding NFPs (13/19) have zero in-house modelling
staff, corroborated by the free text ("No resources for modelling in Luxembourg"; Ireland "limited
capacity"; Italy outsources to Fondazione Bruno Kessler). *Supply:* the RespiCast hubs deliver exactly
what those NFPs lack — a multi-model, ensemble-backed forecast almost every week across two seasons,
covering 15–30 countries. Yet **awareness** of RespiCast stays modest (≈2.1/5, a recency effect) and flat
across capacity bands (and the survey measures awareness, not use). **The binding constraint is on the
demand side** — awareness, interpretation
support, and a route into the decision — not the availability of a forecast.

## The COVID-19 hospitalisation continuity claim

Claim under test: *"COVID-19 hospitalisation forecasts have been produced almost uninterrupted since
2021, first by the COVID-19 hub and then within RespiCast."*

**Verdict: supported, with one refinement.** Stitching the archived EU COVID-19 Forecast Hub to
RespiCast-Covid19 on the ISO-week grid, hospitalisation forecasts ran **2021-07-26 → 2026-06-22**:
**255 of 257 weeks (99.2%)**, only **two isolated one-week gaps**, and a **seamless** handover (the
archive's last hospitalisation round is the ISO week immediately before RespiCast's first). The one
refinement worth stating plainly: **"since 2021" means late July 2021** for hospitalisations
specifically — COVID *cases and deaths* go back to February 2021, but hospitalisations were added mid-
2021. Two honest caveats on the word "produced": (i) in the EU-COVID-hub era the record was **thin** —
a median of ~4 models per week — before roughly doubling under RespiCast; (ii) coverage is *presence*,
not skill. A defensible restatement: *"An ensemble COVID-19 hospitalisation forecast has been published
for the EU/EEA almost every week since mid-2021, without a meaningful interruption across the change of
host from the EU COVID-19 Forecast Hub to RespiCast."*

## Independent verification pass

A separate multi-agent pass re-derived the headline numbers from the **raw** files, blind to
`output/`. Results:

- **Survey — all matched** except one rounding nit: Q10 agree+strongly = 9/19 = 47.4% (reported as
  **47%**, not 48%). The pipeline reports 47%.
- **Hub targets, date ranges, round counts (88 / 91), and country counts (30 / 25 / 15) — all matched.**
- **Per-week model counts — the pipeline matched the raw recompute exactly** (e.g. ILI 2026-03-04 = 16;
  ARI 2026-03-04 = 10; ILI 2024-10-23 = 15). Earlier hand-guessed figures in the verification prompt
  were the ones that were wrong; the pipeline was right.
- **Ever-active model counts (ILI 25, ARI 18, COVID 17) — all matched.**

The one substantive correction the pass forced was narrative, not numeric: an earlier draft claimed
"engagement rises where capacity is thin". The data do not support it (short-term forecasting engagement
is flat: 2.42 / 2.25 / 3.0). The framing was corrected to "near-universal need, flat uptake →
demand-side bottleneck". See `output/audit.json` for the full audit.

## Caveats — read every claim through these

- **n = 19** of ~30 EU/EEA NFPs. Every share moves ~5 points per respondent; the `>10 staff` band is a
  **single** institute. No confidence intervals or significance are reported, by design.
- **Q5 measures AWARENESS, not use** (*"How aware are you…"*, 0–5). Do not translate scores into use,
  uptake, or "% who find X useful". (Internal `eng_*` key names are historical — see `decisions.md`.)
- **The Q5 awareness gradient tracks time in circulation** (guidance 2020, RespiCast 2023, RespiCompass
  2024): newer = less known, not less valuable.
- **Q8 is prospective** ("will be useful") — stated expectation, not demonstrated use.
- **Coverage counts presence, not skill.** Nothing here shows the ensemble beat the baseline or was
  well-calibrated.
- **NFP awareness ≠ national involvement.** Several countries have in-house or partner teams already
  submitting to the hubs (Institut Pasteur, UHasselt/UNamur, Bruno Kessler), so "NFP not aware" can
  understate country-level involvement.
- **Country coverage ≠ national uptake.** A country appearing in a panel means someone forecast it, not
  that its NFP uses the output.
- **Non-response is unaddressed** — the 19 respondents may skew more aware and engaged than the ~11
  non-respondents.

## Highest-value next analyses (ranked)

1. **Bring in forecast skill.** Compute WIS / interval coverage / calibration of the hub ensemble vs
   the quantile baseline and individual models. This is the single biggest missing piece — it converts
   "coverage" into "value", and NFPs explicitly asked for transparent performance metrics.
2. **Directly test the capacity → supply link.** Map each of the 19 NFP countries to whether a team
   from that country submits to the hubs, then cross-tab in-house staff and Q8 value against national
   hub participation — turning two parallel datasets into one joined story.
3. **Build an awareness → value funnel** (aware → engaged → intends to use → has an integration
   mechanism) to locate the demand-side bottleneck quantitatively.
4. **Recency-adjust awareness:** compare RespiCast against the COVID Forecast Hub at equal maturity
   (months since launch) to separate time-in-market from intrinsic value.
5. **Per-country coverage depth:** how many countries have ≥3 models every week vs sporadic coverage?
   Flag thin-coverage countries (COVID reaches only 15) as capacity-building targets.
6. **Characterise non-response** — which countries did not reply, and do they skew high or low capacity?

## Recommendations for ECDC (from the theme synthesis)

Position short-term forecasts (RespiCast) as the primary decision-relevant product and scenarios
(RespiCompass) as explicitly-caveated planning tools; publish transparent skill metrics and consider a
smaller, curated ensemble; provide plain-language outputs and short training for non-modellers
(surveillance specialists, policymakers); integrate forecasts with surveillance (an explicit
RespiCast↔ERVISS link); and back an ECDC recommendation plus sustainable EU funding for national
modelling capacity, so uptake can happen equitably. Full list in `output/themes.json`.

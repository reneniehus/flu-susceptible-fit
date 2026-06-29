# Inter-seasonal influenza variability in the EU/EEA: what phenomenological descriptors of ILI+ curves tell us

*Analysis of 166 country-seasons (25 EU/EEA countries, 8 seasons) using the descriptive (smoothed-curve)
method. Produced by a multi-agent analysis with adversarial verification of every claim; figures in
`output/analysis_patterns.png`, stats in `code/05_analysis/analyse_patterns.R`.*

## Framing

We characterise 166 country-seasons by descriptors read off centred-moving-average smooths of weekly
ILI+ curves: AUC and peak height (burden), onset/peak week (timing), and steepness (exponential growth
of the rising limb). Two distinctions govern everything below. First, **scale**: AUC and peak height
vary ~1000x across countries by reporting scale and are interpretable only *within* country (log +
country-demeaning); timing and steepness are scale-free and comparable across countries. Second, the
**source/era confound**: pre-COVID == RespiCompass, post-COVID == ERVISS, so any era difference is
inseparable from a measurement change except where a single bridging season breaks it.

## 1. Patterns among descriptors (best-supported first)

**Take-off speed is not a strong driver of total burden.** Within country, steepness is uncorrelated
with AUC: pooled r = 0.04 (95% CI crosses zero); the confound-purged estimate (demeaning by
country×source) is ~0.12 and still not distinguishable from zero. AUC is instead dominated by peak
height (within-country r = 0.94) and modulated by duration (AUC~n_weeks = -0.21). The honest framing
is "weak, not distinguishable from zero," **not** "decoupled/zero" — the pooled 0.04 partly averages an
opposite-signed pre (+0.15) and post (-0.17) signal. The better-established companion result is a
**modest, robust peakedness effect**: steepness~peak-height purges to r = 0.25 (p = 0.001, robust to
outlier removal and Spearman). Faster epidemics pile incidence into a sharper, slightly taller peak
without adding total area — and this link is partly definitional, since growth rate and peak are read
off the same rising limb.

**Onset and peak co-move; timing is largely decoupled from intensity.** Within country, onset~peak-week
r ~ 0.58 (survives both eras and joint country×source demeaning). A late season is "sharper in time" —
but **not predictably bigger or smaller**: within-country peak-week~log-peak ~ -0.05 (Spearman ~0), and
any late-onset->intensity slope is trivial. The peak is substantially *anchored*: regressing peak on
onset gives a slope of only ~0.40, so a week-later onset shifts the peak ~0.4 wk and the rest becomes a
shorter lead-in. We **retract** the earlier "rigid calendar translation" and "lateness compresses the
rising limb (r = -0.73)" claims: rise = peak - onset by construction forces that correlation, the
observed value sits *below* an independence shuffle null, and trimming threshold-artifact rows collapses
it to ~-0.27. Honest statement: **anchored peak, variable lead-in.**

**The onset~steepness sign flip is a genuine Simpson reversal — but descriptive, not causal.** Within
country onset~steepness ~ 0 (-0.04); across the 25 country means it is +0.56 (p = 0.004, bootstrap- and
LOO-stable, robust to era and to onset-threshold artifacts). The law of total covariance splits it
cleanly (positive between-country, ~null within), and since ~75% of variance is within-country the
pooled estimate (+0.10) hides it. **Crucially, the across-country +0.56 is near-tautological**: peak
week is nearly fixed across countries, so a later onset mechanically forces a shorter, steeper rise;
partialling out rise time collapses it to -0.10. It is real on average but season-heterogeneous (per-
season -0.13 to +0.66) and should not be read as a country-level transmission law.

## 2. Epidemiological / biological meaning (interpretation, flagged as such)

Mechanistically, growth rate r ~ (R0-1)/generation-interval reflects the take-off *engine* —
transmissibility x starting susceptibility x seasonal forcing — whereas total burden is a saturating
final-size integral plus duration plus reporting scale. These share only R0/susceptibility, which
rationalises the near-zero steepness~AUC link. **A central caveat runs through all timing/steepness
results**: national ILI+ curves are spatial overlays of many local epidemics (Viboud 2006; Gog 2014),
so "steep" national rises partly index spatial synchrony, not local R. The across-country (but *not*
within-country) steepness~onset structure fits this synchrony reading. None of the final-size,
susceptibility-rebuild, or antigenic-match mechanisms can be tested here (no attack-rate, subtype, or
latitude data), so they remain interpretation.

**Post-COVID timing shift.** Onset is earlier post-COVID (20.5->17.2 wk). This *partly* survives the
source confound: the 2023/24 season was measured mostly by RespiCompass, and holding source constant
still shows earlier onset (within-country paired, LOO-stable, p < 0.001), concordant with an internal
ERVISS trend. After removing LU/MT data-quality outliers and properly accounting for the surveillance-
window detection artifact (correctly via country fixed effects, not the pooled model), the defensible
magnitude is **~1.5-2.5 wk earlier onset and ~1.5-2 wk earlier peak** — not 3.9 wk. The **"less steep"
claim is dropped** (within-source NS; ERVISS internally trends steeper). Only one season breaks the
confound, and subtype composition is uncontrolled, so give this moderate weight.

## 3. Vaccination linkage (65+ coverage)

There is **no detectable, robust link** between 65+ coverage and all-age ILI+ descriptors. Within
country every correlation is ~0 (|r| <= 0.08, all p > 0.45); across-country peak-week (-0.33) and onset
(-0.29) are non-significant, below the n=22 detection floor (|r| > 0.42), single-country-driven, and
sign-unstable across seasons. They are also **one weak timing signal counted twice** (onset and peak-
week are collinear; rise length ~ coverage r = -0.01). This is "no *detectable* effect," not "proven
null": ~97% of coverage variance is between-country, so the only contrast with real predictor variation
is the cross-country one — exactly the one confounded by latitude, health system, and age structure.
Interpretively, 65+ are largely on the receiving end of transmission while the all-age curve's
growth/timing is set by higher-contact younger groups, and timing is partly synchrony-driven — neither
of which 65+ coverage can move. (Cleanly, this analysis is 100% pre-COVID/RespiCompass, so it is free
of the source confound.)

## Implications for the broader project

As candidate "sensors," **timing descriptors (onset, peak week) and the steepness~peak-height peakedness
channel are the most promising** — they are scale-free, internally consistent, and carry a genuine
earlier-onset post-COVID signal worth tracking. **Steepness should be treated as a synchrony/timing
quantity, not a burden or local-R proxy**, and any descriptor combining onset, peak week, and steepness
must guard against the rise = peak - onset definitional coupling. Highest-value follow-ups: bring in
subtype (H3N2/H1N1/B), latitude/connectivity, and a sub-national (spatial-synchrony) decomposition to
convert these phenomenological descriptors into mechanistically interpretable sensors, and validate the
post-COVID timing shift as ERVISS accumulates more seasons under a single source.

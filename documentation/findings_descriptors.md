# Inter-seasonal influenza variability in the EU/EEA: what phenomenological descriptors of ILI+ curves tell us

*Analysis of 166 country-seasons (25 EU/EEA countries, 8 seasons) using the descriptive (smoothed-curve)
method. Produced by a multi-agent analysis with adversarial verification of every claim. Figures in
`output/` (`analysis_patterns.png`, `hierarchical_effects.png`, `bayes_subtype.png`,
`bayes_prior_burden.png`); stats in `code/05_analysis/` (pooled correlations in `analyse_patterns.R`,
partial pooling in `hierarchical_models.R`, subtype in `bayes_subtype.R`, prior burden in
`bayes_prior_burden.R`).*

## Within-country hierarchical models (partial pooling) — the precise answers

*Cross-country correlations are dropped as confounded (latitude, health system, age structure). This
section is WITHIN country with partial pooling (lme4: predictor group-mean-centred within country,
country random intercepts + random slopes), so the shared EU/EEA slope learns from all 25 countries
while across-country heterogeneity is estimated. Slopes are standardised (log AUC/peak) ≈ within-
country correlations; 95% CI in brackets; figure `output/hierarchical_effects.png`. Frequentist partial
pooling; a full Bayesian fit would give near-identical inference at this n.*

**Which 2 of {steepness, peak height, AUC} share the most information? → AUC and peak height,
overwhelmingly.** Within-country standardised slope peak-height~AUC = **0.98 [0.90, 1.05]** (near-
collinear, and largely mechanical — peak height is the dominant component of AUC). Steepness shares
almost nothing with burden: steepness~peak-height = 0.07 [0.01, 0.14] (tiny but non-zero, the
"peakedness" channel), steepness~AUC = 0.02 [−0.05, 0.08] (null).

**What does onset week predict? → the peak's TIMING, and nothing about its size.** onset→peak-week =
**0.99 [0.69, 1.28]** (a season starting a week later peaks ~a week later, almost 1:1; real across-
country heterogeneity, random-slope SD 0.61). But onset does NOT predict intensity: onset→peak-height =
0.07 [−0.05, 0.18], onset→AUC = 0.06 [−0.04, 0.16] (both null), nor steepness (−0.05 [−0.33, 0.22]). So
onset is a clean TIMING sensor, uninformative about how big or how steep the season will be.

**What does 65+ vaccination coverage predict? → nothing detectable.** Every within-country effect is
null with CIs spanning zero (pre-COVID block, n=97, 22 countries): peak-height 0.13 [−0.29, 0.55], AUC
0.14 [−0.26, 0.53], onset −0.23 [−1.18, 0.72], peak-week −0.24 [−1.35, 0.87], steepness 0.26 [−0.73,
1.25]. Wide intervals reflect limited within-country coverage variation → "no detectable effect," not
"proven zero." (Random slopes were not identifiable here → random-intercept fits.)

*(The cross-country analysis below is retained for context but is NOT used for inference — confounded.)*

## Dominant subtype and the subtype → descriptor analysis

**Method (WHO/ECDC).** Dominant (sub)type per country-season = the (sub)type with the largest share of
**characterised** influenza detections over the season, among A(H1N1)pdm09, A(H3N2) and B (sentinel +
non-sentinel combined), requiring ≥ 20 typed detections for a call. 84% of country-seasons have a clear
plurality (top share ≥ 50%; median 64%).

**Coverage.** ERVISS typing starts 2021, so only the post-COVID analysis seasons get a subtype: **59 of
166 descriptor country-seasons (36%)** — all of 2023/24, 2024/25, 2025/26; none of the pre-COVID seasons
(2014/15–2018/19), which would need external WHO/ECDC season reports.

**Critical confound.** The dominant subtype is essentially **season-determined**: 2023/24 = H1N1 (23/29
countries), 2024/25 = B (23/31), 2025/26 = H3N2 (24/30). With only three post-COVID seasons, subtype is
nearly collinear with season, so a "subtype effect" cannot be separated from anything else specific to
those three seasons (vaccine match, prior-season immunity, weather). Read the results as among-**season**
differences labelled by subtype, NOT causal subtype effects.

**Bayesian model** (Gibbs sampler: descriptor ~ subtype + country random intercept, weak priors;
posterior means match lme4; R-hat ≤ 1.01). Within-country contrasts (SD units; * = 95% CrI excludes 0;
figure `output/bayes_subtype.png`):
- **Burden:** B-dominant seasons have higher AUC than both H1N1 (+0.31*) and H3N2 (+0.22*), and higher
  peak height than H1N1 (+0.24*).
- **Timing:** H3N2 seasons peak ~1.1 SD **earlier** than H1N1 (−1.11*); B peaks and onsets later than
  H3N2 (peak +1.36*, onset +0.78*). Ordering H3N2 → H1N1 → B from earliest to latest.
- **Steepness:** no subtype difference (all CrIs span 0).
These directions align with known epidemiology (H3N2 seasons earlier/sharper; B later) but, given the
confound, remain among-season differences pending more single-source seasons.

## Prior-season burden (immunity carryover)

Does the previous season's burden predict the current season, within country? Burden is lagged by
exactly one year (the COVID gap auto-excludes invalid lags) → **120 lagged country-seasons** (25
countries; 2015/16–2018/19 and 2024/25–2025/26). Bayesian hierarchical model borrowing strength across
all countries (Gibbs: country random intercept **and** random slope; shared EU/EEA slope; R-hat ≈ 1.00,
matches lme4). Predictor and outcomes standardised, predictor within-country-centred.

**Prior burden is a weak-to-null predictor of everything** (all 95% CrIs span zero, SD units): current
AUC −0.06 [−0.39, 0.24], peak height −0.08 [−0.41, 0.21], peak week +0.36 [−0.41, 1.05], onset +0.13
[−0.57, 0.77], steepness +0.03 [−0.61, 0.68]. The point estimates lean *faintly* immunity-consistent
(bigger prior → slightly smaller burden, later peak), and notably the AUC/peak slopes are **negative
despite** a reporting-persistence confound that should push them positive — so a small true depletion
signal is plausible but not credible at this n. The prior→peak-week effect is also very heterogeneous
across countries (random-slope SD 0.79).

Why so weak — and the proper test we can't yet run: (1) **subtype rotation breaks cross-protection** —
total prior burden is a poor immunity proxy when the dominant subtype changes year to year (it does:
H1N1→B→H3N2 post-COVID), so the right predictor is prior burden *of the same subtype*, which our data
can't support yet (subtype is post-COVID only, and those seasons rotate); (2) antigenic drift restores
susceptibility annually; (3) reporting persistence masks the burden magnitude; (4) national curves are
spatial overlays. Figure `output/bayes_prior_burden.png`.

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

## Pre-COVID VE + dominant subtype (first single-era model)

`code/05_analysis/bayes_precovid_ve_subtype.R`. **Pre-COVID only** (2014/15-2018/19: 107 country-seasons,
23 countries, 5 seasons; 97 with 65+ coverage) — one clean era, avoiding the pre/post-COVID
behaviour/reporting confound (see `reflections.md`). Within-country partial pooling (country random
intercept), Gibbs, lme4-validated (R-hat approx 1.00). Outcomes standardised: log(AUC), onset week,
log(peak height). **Model 1** (as specified): dominant subtype (categorical) + VE (season-level VE against
the dominant subtype). **Model 2** (mechanistic refinement): subtype + **protection = VE x 65+ coverage**,
which varies by country x season. Effects in SD units, 95% CrI; `*` = excludes 0.

**Robust signal — subtype contrasts (consistent with the 8-season analysis):** **B seasons carry the
largest burden** (AUC B-H1N1 **+0.50***, B-H3N2 **+0.63***) and highest peak; **A(H3N2) onsets earliest**
(H3N2-H1N1 **-1.54*** on onset week), B latest (B-H3N2 **+1.86***). *Caveat:* pre-COVID **B rests on a
single season** (2017/18), so its contrasts are essentially "2017/18 vs the rest" — season-confounded.

**VE / protection — the predicted vaccine mechanism is NOT evident, and a confounded one appears.** The
mechanistic prediction (see `reflections.md`) was: protection should push **burden (AUC/peak) down** and
leave **timing (onset) unmoved**. Instead:
- burden: protection slope on AUC **-0.08 (ns, -0.25 to 0.09)**, on peak **-0.03 (ns)** — the expected
  reduction is weak and not distinguishable from zero;
- timing: the **only** significant VE/protection association is with **earlier onset** (VE **-0.65***;
  protection **-0.38***) — but vaccines cannot plausibly shift onset timing, so this is **season-level
  confounding**: VE co-varies with subtype and with season-specific onset. Even the country x season
  protection term is dominated by its season-level VE component (65+ coverage is nearly constant within a
  country), so it **inherits the confounded timing signal instead of isolating the burden mechanism**.

**Takeaway (the identifiability limit made explicit):** with subtype and VE both **season-level** over only
5 seasons, and coverage carrying too little within-country variation to rescue it, the vaccine-**burden**
effect is not identifiable here and apparent "VE effects" are entangled with subtype/season. The subtype
contrasts are the informative — if season-confounded — signal. This is the clean pre-COVID baseline against
which a later, era-controlled model that adds the post-COVID seasons should be read.

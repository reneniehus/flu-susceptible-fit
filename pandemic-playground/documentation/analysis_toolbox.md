# The analysis tool box — phases, methods, truth, extension fronts

The tool box is organised by the **phases of a response**. Each tool consumes the tidy `observed`
schema (`as_analysis_input()`, or real data in the same shape), returns an estimate, and has a
`score_*()` companion that grades it against the matching quantity in `truth`. Run
`Rscript demo/run_playground.R` to see all ten scored at once.

A recurring theme: several tools are **exactly correct on the latent truth** but **biased on the
observed data** — and the bias is the lesson. The playground exists to show these, not to hide them.

**Not just COVID.** Where a tool has a timescale (the Cori window, the growth window, the forecast
horizon, the SIRS recovery rate), its default is *derived from the generation interval* rather than a
fixed COVID-era calendar constant, so the tools track the pathogen — fast flu to slow measles. And where
signal is too thin to support a method (few cases, a burnt-out tail, a growth-rate CI straddling zero),
the tool **refuses** (NA / empty) rather than emit a prior-dominated confident-but-wrong number.

---

## Phase 0 — before local introduction (the source X still dominates)

### Is it spreading, and how fast? — `phase0_growth_R`
Poisson-GLM growth rate `r` (doubling time `ln2/r`), converted to `R` through the generation interval
by Euler–Lotka / Wallinga–Lipsitch, `R = 1 / Σ_a g(a) e^{−r a}`.
**Scored against** the true early `R`. **Known bias:** on the *observed* case curve during scaling-up
testing, `r` (and hence `R`) is inflated by the ascertainment ramp; the same tool is exact on the true
infection curve. **Extension:** `EpiEstim` for R directly (the Cori tool below applies to the same
early curve).

### Is it spreading, and how *heterogeneously*? — `phase0_clusters`
The one early method that recovers `R` **and** the superspreading dispersion `k` together, from the sizes
of observed transmission chains. A cluster of final size `n` under negative-binomial offspring has the
Blumberg–Lloyd-Smith chain-size probability (via the Borel/Otter identity), so a joint MLE on the observed
sizes returns `(R, k)` with CIs; a Poisson-limit `k→∞` recovers the classic Borel distribution.
`cluster_size_fit` (the estimator), `cluster_extinction` (chance a chain dies out), `cluster_singleton_k`
(the singleton-fraction shortcut), `cluster_analysis` (end to end on `observed$clusters`).
**Scored against** the true `(R, k)` used to grow the chains. **Known subtlety:** a subcritical `R<1` and a
supercritical `R>1` can imply *near-identical* cluster-size distributions (the Nishiura conjugate
ambiguity), and `k` governs how sharply the sizes pin `R` at all — so few, small clusters leave `R` genuinely
uncertain, which the tool reports rather than hides. **Extension:** `{epichains}` for the exact likelihood
and simulation; phylodynamics (`{BEAST}` skyline) for `R`/`k` from genomes.

### How big is the source really — and how much are we missing? — `phase0_catchment`
Back-calculate the source's infectious prevalence from exported cases and travel fractions,
`prevalence_X ≈ detected exports / Σ_c (travel-fraction_c · detection_c)`, anchored on well-surveilled
destinations — the method that sized Wuhan (~40× the reported count; Imperial / Wu-Leung-Leung) and Mexico
2009 H1N1 (Fraser et al.) from cases detected *abroad*, independent of the source's own thin surveillance.
`catchment_backcalc` gives the point estimate; `catchment_range` returns the honest **low/central/high**
range across two real-world corrections — visitor-vs-resident representativeness `V = 1 − e^{−(r+γ)d}`
(De Salazar; short-stay travellers under-sample a growing source) and differential detection (Niehus;
benchmark anchors against the best surveillance). Both corrections only ever push the estimate **up**.
**Scored against** the true source prevalence (and the true source under-ascertainment). **Extension:**
a small Poisson model `detected ~ Poisson(volume · detection · prevalence)` in place of the arithmetic.

### Who is under-detecting their imports? — `phase0_importation_risk`
Regress detected imports on flight volume across the well-surveilled anchor countries (De Salazar),
predict expected imports for all, flag those below the prediction interval.
**Scored against** true surveillance quality (flagged countries should be the poorly surveilled ones).
**Extension:** negative-binomial regression if the imports are overdispersed.

---

## Phase 1 — early local exponential growth

### How deadly is it? — `phase1_cfr`
Delay-adjusted (confirmed) CFR (Nishiura): `cCFR(t) = deaths(t) / Σ_o cases_o · F_death(t−o)`;
`cfr_static` / `cfr_rolling`.
**Scored against** the true confirmed CFR (eventual deaths ÷ eventual cases), and reported next to the
true IFR. **Known bias:** the adjustment removes the *delay* bias but not case *under-ascertainment*, so
the cCFR overstates the IFR by ~`death_detection / ascertainment` — the gap the playground displays.
**Extension:** `{cfr}` (`cfr_static`/`cfr_rolling`); IFR from seroprevalence once available.

### Growing or shrinking, right now? — `phase1_rt`
Self-contained Cori/EpiEstim renewal Rt: conjugate `Gamma(a+ΣI, 1/b+ΣΛ)` over a sliding window.
**Scored against** the realized instantaneous `Rt_effective`. **Known biases:** feeding onset-dated
cases lags infection-time R by ~one incubation period; importation contaminates a country's early R
(imports read as local); intervals under-cover a stepped truth near change points. Exact on the
import-free source's true infections. **Extension:** `EpiNow2` (infers infections; delay- and
truncation-aware).

### How many recent cases, really? — `phase1_nowcast`
Invert the reporting-delay truncation on the triangle, `nowcast_o = observed_o / F(as_of − o)`, with an
assumed or empirically-estimated delay and a honesty guard that flags days too recent to nowcast.
**Scored against** the eventual onset totals (vs just trusting the truncated data). **Extension:**
`epinowcast` / `EpiNow2::estimate_truncation` for the full Bayesian version.

---

## Phase 2 — established transmission, growth to peak, healthcare demand

### Will we breach capacity? — `phase2_forecast`
Nowcast the recent case curve, estimate current R (Cori), project cases forward by a fixed-R renewal
over an R-scenario set, map to admissions via a fitted case→admission ratio and the onset→admission
delay, compare to capacity.
**Scored against** the admissions that actually occurred (breach call and timing). **Known bias:** a
fixed-R forecast cannot see a *new* change in transmission and overshoots at the turnover; the onset-R
lag delays detection of the peak. Good in the clear growth and decline phases. **Out-of-envelope flag:**
because there is no behavioural feedback, the projection can shoot straight *through* capacity — which real
systems don't do — so a scenario whose peak exceeds `implausible_multiple × capacity` is flagged
`out_of_envelope` and read as "capacity at serious risk", never as a literal number (see the well-posedness
note below). **Extension:** `EpiNow2` forecast; a time-series foundation model / method-of-analogues where
the mechanistic overshoot is least trustworthy.

### Are controls working? — `phase2_intervention`
Interrupted time series: a segmented Poisson model with the growth rate breaking at the intervention
date (aligned by the onset lag), converted to a before/after R.
**Scored against** the realized R either side of the intervention. **Reported as association, not
proof** — seasonality, behaviour and depletion can bend the curve simultaneously. **Extension:**
`EpiNow2` breakpoints / `EpiEstim` windows.

---

## Phase 3 — sustained transmission, later waves, endemic-ish

### Will a variant take over, and how fast? — `phase3_variant_selection`
Binomial GLM of the sequenced counts, `logit P(variant) = a + s·day`; `s` is the selection coefficient,
with the 50% crossover day.
**Scored against** the realized logit-slope of the true variant frequency. **Known bias:** with large
sequenced counts the CI is very tight and can under-cover, since the true `s` drifts as R changes.
**Extension:** `nnet::multinom` for >2 co-circulating variants.

### Next season, and does boosting help? — `phase3_scenarios`
A parsimonious SIRS integrated with base-R RK4, run over a waning × transmissibility × booster-uptake
factorial; outcomes (peak, timing, cumulative incidence) reported relative to a reference and ensembled across
the waning assumption; seed the initial immune fraction from where the current run left the population.
**Validated** against the analytic SIR final-size relation. This is scenario *exploration*, not point
estimation, so it is checked by the integrator's correctness rather than scored against a single truth.
**Least well-posed tool in the box, deliberately:** over a season, behaviour, variant properties and waning
are assumed not measured, so use it to *rank levers* under shared assumptions, never to quote an absolute
next-season number. **Extension:** `deSolve` / `odin` SIRS, optionally age-structured; a Scenario-Hub target format.

---

## Which questions are well-posed — and where foundation models come in

Not every question a compartmental model *can* be asked is one it can answer honestly. It is worth being
explicit about the gradient, because the tools above span it:

- **Well-posed (short horizon, inside the observed envelope).** Logistic variant selection (`s` is
  identifiable from sequence counts and extrapolates cleanly for weeks), reporting-triangle nowcasting
  (a bounded deconvolution of a delay we can estimate), and near-term Rt. These recover a quantity that is
  actually *in* the data. Score them, trust them, correct their biases.
- **Fragile (medium horizon).** The capacity forecast. Well-posed for one to two weeks of clear
  growth/decline; **misleading** the moment it implies admissions shooting *through* a threshold, because
  the mechanistic model has no behavioural-feedback term and real hospitalisations don't do that — hence the
  `out_of_envelope` flag. Read it as a risk signal past that point, not a number.
- **Ill-posed as prediction (long horizon).** Next-season booster scenarios. The mechanistic detail is
  swamped by assumptions (behaviour, the next variant, waning) no one can measure in advance. Legitimate
  only as *relative* lever-ranking, never as an absolute projection.

**Foundation models are the natural counterpart — not a blanket antidote.** A mechanistic model is
overconfident *outside* its envelope (it happily extrapolates its own equations). A data-driven / time-series
foundation model — method of analogues, seasonal-MoA, TabPFN-TS, PandemicLLM, and mechanistic-statistical
hybrids — is instead grounded in *how real curves have actually bent*, so it is more honest exactly where the
compartmental forecast overshoots. But it is symmetrically overconfident *inside* its training distribution
and untrustworthy when the new pathogen is unlike anything it has seen. "Antidote" is the wrong metaphor:
the two fail in opposite regimes. The real prize is **honest out-of-distribution detection** — a model that
says *this is unlike what I was built on* — which is why the tools here carry explicit envelope flags rather
than pretending the projection is always a prediction. (Fuller argument in `reflections.md`.)

---

## The extension fronts in one place

| Front | Where | What drops in |
|---|---|---|
| Curated delay distributions | `epidist.R::as_epidist()` | a real `{epiparameter}` object |
| Cluster / branching R and k | `phase0_clusters` | `{epichains}` (exact likelihood + simulation); phylodynamic skyline |
| R estimation | `phase1_rt`, `phase0_growth_R` | `EpiEstim`, then `EpiNow2` |
| CFR | `phase1_cfr` | `{cfr}` |
| Nowcasting | `phase1_nowcast` | `epinowcast`, `EpiNow2::estimate_truncation` |
| Forecasting | `phase2_forecast` | `EpiNow2`; time-series foundation models / method-of-analogues where the mechanistic overshoot is least trustworthy |
| Variant models | `phase3_variant_selection` | `nnet::multinom` for >2 variants (not implemented), phylodynamics |
| Scenario ODEs | `phase3_scenarios` | `deSolve` / `odin`; age structure |
| Age-structured DGP | `renewal.R`, `draw_parameters.R` | age classes with a contact matrix |
| Partial cross-immunity | `renewal.R` | strain-specific susceptible pools |
| Out-of-distribution detection | all forecast/scenario tools | an explicit "is this input unlike what the method was built on?" check |
| Real data | `analysis_common.R` | real surveillance frames in the schema (see `real_data.md`) |

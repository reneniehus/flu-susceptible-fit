# Design decisions & rationale

The README records *what* the playground does; this file records *why* the key design choices were
made, so the reasoning survives beyond commit messages. Each entry: the **decision**, the reason, and
the main alternative considered. Append new decisions as they are made.

## Architecture

- **Strictly separate the latent data-generating process from the observation model, and return
  BOTH.** This is the whole reason the playground exists: an analytical tool can only be *scored* if the
  truth it is trying to recover is available alongside the degraded data it is given. So the pipeline is
  two clean halves — the DGP (`simulate_source/flights/importations/local`) produces `truth`, and
  `observe()` degrades it into `observed` — and no analysis tool is ever allowed to read `truth` except
  through a `score_*()` function. *Alternative:* generate only observed data (a plausible simulator) —
  rejected, because then nothing can be scored and the instrument loses its point.

- **One config object plus one seed fully determine a run.** Every assumption lives in
  `default_config()`; `set.seed(config$seed)` is called once, in `simulate_pandemic()`. An experiment is
  an edit to a copy of the config. *Alternative:* scattered constants / multiple seeds — rejected as
  irreproducible and hard to reason about.

- **The config is a plain nested list, not an S4/R6 object.** Trivial to read, print, edit in a script
  and serialise; matches the parent repo's `settings()`. `validate_config()` gives the safety an object
  system would, without the ceremony.

## The delay distributions and the emulated package stack

- **Delay distributions are a tiny, {epiparameter}-SHAPED local class (`epidist`), not the
  {epiparameter} package.** The brief asks for delays "as {epiparameter} objects". {epiparameter} is the
  right home for curated, citable distributions, but it is not installable here (CRAN is blocked in the
  managed environment). So `epidist` mirrors its shape — a name, a family, the family's parameters — runs
  on base R, and **accepts a real `<epiparameter>` object via `as_epidist()`** the moment the package is
  present. Everything downstream calls only `discretise()`, so the swap is a one-line change at the
  config. *Alternative:* depend on {epiparameter} — rejected (unavailable; also adds a heavy dependency
  to an otherwise dependency-free engine).

- **Emulate the whole production tool stack (EpiEstim, {cfr}, EpiNow2, epinowcast, deSolve/odin, nnet)
  leanly in base R, rather than depend on it.** Three reasons: CRAN is blocked (only the tidyverse,
  `nnet`, `MASS`, `mgcv` are present); a self-contained engine is more reproducible and portable; and
  implementing each method transparently makes its **assumptions and biases visible** — which is exactly
  what a testing playground should surface, not hide behind a package call. Each tool names the package
  it emulates and the extension front where the fuller version drops in. *Alternative:* wait for the
  packages — rejected (blocks all progress and hides the mechanics). This mirrors the parent repo's own
  choice to write its SIR/EKF rather than import a package.

- **Discretise continuous delays by interval censoring, with two day-0 conventions.**
  `P(d) = F(d+½) − F(d−½)`, renormalised (Cori 2013). The generation interval uses the `cori` boundary
  (`P(0)=0`: no same-day transmission, matching EpiEstim's `w_0=0`); incubation / reporting / onset→death
  keep the day-0 mass (`interval`), where a zero-day delay is real. *Alternative:* a single convention
  for all — rejected as wrong for the generation interval.

## The data-generating process

- **Discrete-time daily renewal WITH susceptible depletion, one engine for the source and every
  country.** `I_t ~ Poisson(R_t · S_{t-1}/N · Σ_s g(s) I_{t-s})`. The depletion factor `S/N` is essential:
  without it epidemics never peak, and Phases 2–3 (growth to peak, healthcare demand, later waves) have
  nothing to work on. Using the *same* engine for source and destinations means any tool that works on
  one works on the other, and the source and its imports are generated identically. *Alternative:* a
  pure exponential/branching process (no depletion) — rejected (no turnover); a full ODE SIR — heavier
  and awkward to seed stochastically from a handful of imports.

- **The generation interval doubles as the infectiousness profile; infectious PREVALENCE is incidence
  convolved with its survival.** In renewal theory `g(·)` *is* the normalised infectiousness profile, so
  the probability of still being infectious `s` days after infection is its survival `1−CDF(s)`, and
  prevalence (the currently-infectious count that drives exports) is `Σ_s I_{t-s}·(1−CDF(s))`. This needs
  no extra parameter. *Alternative:* a separate infectious-period distribution — cleaner in principle,
  but an extra knob the brief does not list; noted as an extension front.

- **The two strains share one susceptible pool (complete cross-immunity).** The simplest coupling that
  still lets a fitter variant compete and win. *Alternative:* partial cross-immunity / immune escape —
  richer and realistic, but more parameters; a documented extension front.

- **The variant travels abroad through imports, not by separate per-country seeding.** Imports from X
  are split into wild-type / variant by a binomial draw on X's current variant *prevalence*, so local
  variant emergence follows from connectivity and the source dynamics with no extra machinery. This
  couples the source variant sweep to every country's variant sweep in one mechanism. *Alternative:*
  seed the variant independently in each country — rejected (arbitrary, and severs the link the whole
  importation story is about).

- **Exports ~ Poisson(true flight volume × source prevalence fraction).** The cleanest expression of
  "travellers carry infection in proportion to how infectious the source currently is". Volumes scale
  with destination population (a power law) times a persistent per-country connectivity, modulated by
  annual seasonality and a configurable travel-ban shock. *Alternative:* a full gravity model — more
  realism than a testbed needs.

## The observation model

- **Degrade at the AGGREGATE level with Poisson / multinomial thinning, not per-infection
  microsimulation.** With tens of millions of infections, simulating each individual's onset, report and
  death is infeasible. Instead: expected onsets = infections ⊛ incubation; reported cases ~ Poisson(onsets
  × ascertainment); the reporting delay scatters each onset-day's cases across report days by a
  `rmultinom` draw (the triangle); deaths ~ Poisson(IFR × infections ⊛ (incubation⊛onset→death) ×
  death-detection). Fast, exact in expectation, and every quantity stays an integer count. *Alternative:*
  agent-level microsimulation — rejected (does not scale; unnecessary for the questions asked).

- **`observed$cases_by_onset` is RIGHT-TRUNCATED; the eventual onset total lives in `truth` as the
  nowcast target.** The analyst-facing onset curve only contains reports that have arrived by `as_of`, so
  recent onset days look artificially low — exactly the truncation a nowcast must invert. The eventual
  (fully-reported) onset total is the truth to score the nowcast against, so it is kept in
  `truth$cases_by_onset`. `cases_by_report` (by report date) is always complete to `as_of`.
  *Alternative:* expose only the eventual onset curve — rejected (hides the truncation the tool exists to
  fix).

- **Deaths (and admissions) are drawn as independent thinnings of infections, not as a strict subset of
  detected cases.** This keeps the observation model a set of clean, independent Bernoulli/Poisson
  thinnings. **Consequence, stated plainly:** the confirmed CFR that a delay-adjusted estimator recovers
  equals `IFR × death_detection / case_ascertainment`, which *overstates* the IFR whenever case
  ascertainment is below the death-detection rate — precisely the real early-pandemic bias, and one the
  playground then makes visible by also reporting the true IFR. *Alternative:* make deaths a subset of
  confirmed cases — arguably more realistic for line-list data, but it entangles the two thinnings and
  hides the ascertainment bias we want to demonstrate; revisit if a line-list mode is added.

- **Case ascertainment RAMPS UP over time (testing scaling), deaths are detected near-completely.** This
  is the single most important observation-model choice for realism: a scaling testing regime inflates
  the *observed* growth rate and Rt (early cases missed), which is why the Phase-0 growth tool over-reads
  R on observed cases even though it is exactly correct on true infections. Deaths, better ascertained
  and constant, anchor the CFR. *Alternative:* constant ascertainment — rejected (removes the headline
  early-pandemic bias the playground should teach).

- **Store the REALIZED instantaneous reproduction number `Rt_effective` in truth, alongside the nominal
  schedule.** A Cori/EpiEstim estimate off the case curve targets the *realized* instantaneous R — the
  strain-mix- and depletion-weighted `Σ_k R_k S/N Λ_k / Σ_k Λ_k` — not the nominal per-strain step
  function (which ignores the variant and depletion). Scoring against the nominal schedule would unfairly
  penalise a correct estimate. *Alternative:* score against the nominal schedule — rejected as the wrong
  target.

## The analysis tool box (and its honest limits)

- **Every tool consumes a fixed tidy SCHEMA via `as_analysis_input()`, never the sim object.** This is
  the real-data seam: supply the same frames from real surveillance data and the identical tools run
  unchanged (tested in `test-real-data.R`). *Alternative:* tools reach into the sim — rejected (couples
  analysis to simulation and forecloses real data).

- **Keep, and pin as tests, the biases the lean methods exhibit** — they are the findings, not bugs:
  - *Ascertainment-ramp inflation* of the Phase-0 growth-based R (correct on true infections; high on
    observed cases).
  - *Importation contamination* of a country's early Cori Rt (imports read as local transmission → biased
    high while imports dominate).
  - *Fixed-R forecast overshoot at the turnover*, compounded by the onset→infection lag (R estimated from
    onset-dated cases lags infection-time R by ~one incubation period).
  - *Over-tight credible intervals* under model misspecification (Cori assumes Poisson; huge counts give
    narrow intervals that under-cover a stepped truth near change points).
  Each is documented in the tool's header and, for the first three, exercised by a test.

- **The Phase-2 forecast drives off CASES (which lead), then maps to admissions via the admission
  delay — it does not forecast admissions directly.** Estimating R from admissions is lagged by the
  admission delay and overshoots at the peak even worse. Cases lead admissions, so the case-driven
  forecast peaks on time and the "baked-in" admissions (from already-infected people) are put back
  explicitly. *Alternative:* renewal directly on admissions — simpler but the timing is wrong; rejected.

- **Report Phase-3 scenario outcomes RELATIVE to a reference and ENSEMBLE over the most uncertain axis
  (waning).** Absolute next-season numbers depend on parameters no one knows in advance; relative
  differences between scenarios sharing those unknowns are robust, and averaging over waning states the
  uncertainty instead of hiding it. Matches the scenario-hub target format. *Alternative:* headline
  absolute projections — rejected as falsely precise.

- **Effective IFR from the age table; the DGP is not age-stratified.** An age-structured IFR is collapsed
  to a single effective rate by the population age weights (exact for a non-age-stratified epidemic).
  Full age-structured transmission is an extension front. *Alternative:* age-stratify the renewal — a
  large addition the brief lists as optional.

## Generalising beyond COVID (driven by a multi-pathogen stress test)

A multi-agent stress test ran the playground as five very different respiratory pathogens (1918-severe
flu, mild-iceberg flu, SARS-1 superspreader, sub-critical MERS, measles-like) and probed the structural
limits. It confirmed the DGP itself is robust across pathogen speeds and severities, but surfaced one
big missing capability and a family of COVID-tuned analysis defaults. The changes:

- **Transmission overdispersion / superspreading (`dispersion_k`).** The single biggest generalisation.
  The renewal draw was strictly Poisson (variance = mean), so k-type offspring heterogeneity — the
  defining feature of SARS (k ~ 0.16), MERS and even COVID — was structurally unrepresentable, and every
  imported chain established far too reliably (verified: one seed at R0 = 2.4 establishes 88% of the time
  under Poisson vs ~20% under NB k = 0.16). We add one config scalar `dispersion_k` (Inf = Poisson,
  backward-compatible) and draw endogenous infections as `rnbinom(size = k * force, mu = mu)`. Scaling
  `size` with the day's infectiousness (the effective number of infectors) makes the daily total equal
  the exact sum of per-infector NB(k) offspring (Lloyd-Smith et al. 2005): heavily overdispersed when few
  are infectious (imported chains fizzle) and back to Poisson at the peak. *Alternative:* a fixed
  `size = k` — rejected (wrong aggregate variance; would not collapse to Poisson at scale).

- **Derive analysis-tool timescales from the generation interval, not the calendar.** Nearly every tool
  default was a COVID-tuned constant that misbehaves on faster/slower pathogens: a 7-day Cori window is
  2.7 generations for flu (over-smoothed) but 0.6 for measles (unstable); a 0:30 growth window spans a
  fast epidemic's whole rise-and-fall. So the Cori/forecast/growth windows now default to `NULL` and
  derive from the GI (`~1.5 x GI mean`), and the demo passes the SIRS `gamma` as `1/GI_mean`. The current
  numbers remain available as explicit overrides. *Alternative:* keep fixed defaults — rejected (silently
  wrong off the COVID timescale, the exact trap a "new pathogen" playground must avoid).

- **Low-signal guards so tools refuse rather than emit prior-dominated nonsense.** With few cases the
  Cori posterior collapses onto its prior mean (~5), giving a confident R ~ 5 on a sub-critical or
  burnt-out epidemic. `estimate_rt_cori` now returns a typed empty frame (never a bare NULL) and skips
  windows below a small incidence floor — matching the growth / variant tools, which already `stop()` on
  too little signal. The growth doubling-time CI now returns NA (not a spurious finite number) when the
  growth-rate CI straddles zero, and reports a halving time for a declining curve.

- **The confirmed CFR may exceed 1, and its interval must stay consistent.** Because deaths and cases are
  thinned from infections independently, the confirmed CFR = IFR·death_detection/ascertainment
  legitimately tops 1 whenever death detection beats case ascertainment (any high-severity or early run).
  The point estimate is no longer clamped, the interval is a Poisson-count ratio CI that can also exceed 1
  (so we never show a 6.2 estimate beside a [1,1] interval), and the tool warns when the CFR > 1 (a real
  diagnostic that ascertainment is below death detection) and returns NA below a minimum resolved-case
  count. *Alternative:* clamp to [0,1] — rejected (hides the very ascertainment signal the playground
  exists to expose).

- **Susceptible depletion in the forecast projection.** A fixed-R renewal with no S/N term projects
  unbounded exponential growth, so a fast/high-R pathogen forecast implied more admissions than the whole
  population (measles: 42M admissions in 28 days). `renewal_project` now takes a susceptible pool
  (population x an optional serology-based susceptible fraction), depletes the projected mean by the
  remaining fraction, and caps cumulative projected incidence at the pool. The case->admission ratio is
  fitted on an adaptive window (settled days that actually contain admissions), so a fast epidemic no
  longer silently forecasts zero admissions. *Alternative:* forecast admissions directly with no cap —
  rejected (physically impossible projections).

- **Flag depletion-driven declines in the intervention analysis.** On a fast pathogen a control measure
  can land after the epidemic has already peaked from susceptible exhaustion; `intervention_its` now flags
  `depletion_suspected` when the intervention is at/after the observed peak or growth was already negative,
  so burnout is not mis-read as the intervention working.

- **`attack_rate` was a misnomer; renamed to `cumulative_incidence`.** The SIRS accumulator counts
  infection *events* and, under the waning the scenarios explore, exceeds 1 through reinfection — so it is
  a burden measure, not a fraction-ever-infected. Renamed and documented; it still equals the final size
  for a pure SIR.

**Not built (documented extension fronts).** A `symptomatic_fraction` to separate a true asymptomatic
iceberg from surveillance under-detection; sub-daily integration (rejected — the daily engine is an exact
discrete renewal, verified accurate down to GI mean ~2 d, not an Euler approximation); GI-uncertainty
propagation into the r->R interval.

## Infrastructure

- **EU/EEA populations are the only "real" geography; everything else is simulated.** Population drives
  flight connectivity and each country's susceptible pool, so it is the one datum the simulator needs; it
  is a small committed table (27 EU + IS/NO/LI), no external data. *Alternative:* pull real flight /
  contact matrices — rejected (external data, egress-blocked here, and unnecessary for a testbed).

- **No external data files, no epidemiological-package dependencies; tests simulate their own data.**
  The engine is base R; the tool box adds only the tidyverse (the `nnet::multinom` >2-variant path is a
  documented extension, not implemented). Tests build their fixture by calling the simulator, so they run
  offline with nothing to download. Mirrors the parent repo's offline, snapshot-driven testing, but here
  the "snapshot" is generated deterministically.

- **Surveillance quality is exposed to the analyst as a KNOWN covariate, and the two importation tools
  are "known-covariate" tools by design.** `observed$detected_imports` carries each country's exact
  `surveillance_quality`. This is deliberate and standard: the De Salazar / Fraser catchment logic
  *requires* the detection probability to be known at the well-surveilled anchor countries, and in a real
  response this is the role played by an external surveillance-capacity index (e.g. IDVI) or a
  seroprevalence-anchored ascertainment estimate — a covariate the analyst genuinely has. So
  `catchment_backcalc` legitimately divides by it, and `score_importation_risk` is best read as a
  *consistency check* (do the flagged countries line up with the known covariate?) rather than an
  out-of-sample truth comparison. **Caveat, stated plainly:** because the covariate is exact rather than a
  noisy proxy, these two tools' scores are flattering relative to the others (which recover a genuinely
  hidden quantity). Degrading the covariate to a noisy proxy — and having the catchment tool *estimate*
  detection — is a documented extension front. *Alternative:* keep detection fully latent — rejected,
  because then the catchment method has no anchor and cannot be posed at all.

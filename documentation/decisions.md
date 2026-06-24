# Design decisions & rationale

The README records *what* the repo does and how to run it; this file records *why* the key
modelling, method and data choices were made, so the reasoning survives beyond commit messages and
chat. Append new decisions as they are made. Keep each entry short: the **decision**, the reason,
and the main alternative considered.

## Inference & modelling

- **Fix R0 from the literature (`susc_R0 = 1.5`); fit a single susceptibility `S0` per season.** How
  fast a season's virus spreads is set by the effective reproduction number `R_eff = R0 * S0`
  (equivalently the early growth rate `r = gamma*(R0*S0 - 1)`): intrinsic transmissibility `R0` and
  the susceptible fraction `S0` are multiplicatively linked, so a single curve identifies their
  *product*, not each separately — fitting both jointly is degenerate (R0 drifts to its prior and S0
  follows). We therefore fix R0 and let one per-season `S0` act as the "sensor" of how easily that
  season's virus spread. This is appropriate because the factors of interest here act mainly on the
  *susceptible fraction* — vaccine coverage, antigenic match/mismatch and prior-season population
  immunity (see `PROJECT_SCOPE.md`, refs 4, 6–8) — whereas the determinants of intrinsic
  transmissibility (contact patterns, behavioural restrictions, viral subtype infectivity) are either
  out of scope or reasonably treated as season-invariant. Inferring a susceptibility state from
  incidence under a fixed transmission model is the established *susceptible-reconstruction* approach
  (Bjørnstad et al. 2002; Finkenstädt & Grenfell 2000). **Caveat:** holding R0 fixed means any genuine
  season-to-season variation in transmissibility — notably subtype infectivity, with A(H3N2) seasons
  more transmissible than A(H1N1) (Infection, Genetics and Evolution 2018) — is absorbed into the
  fitted S0. The fitted S0 is therefore a *composite* index of "how easily the season spread", not a
  pure immunity measure; in particular, an apparent association between S0 and dominant subtype may
  partly reflect this absorption rather than susceptibility alone.

- **Fix the seed `I0` (`susc_seed_i0 = 1e-5`), plant it at the season start (~Aug), and read the
  reporting fraction `c` only in relative terms.** The seed size `I0` and the reporting fraction `c`
  both act on the *scale* of the modelled curve — a larger seed and a larger reporting fraction are
  largely interchangeable in their effect on the observed level — so they are not jointly identifiable
  from one curve. We resolve this by fixing `I0` (a small constant import, encoding the assumption
  that flu is reliably re-seeded each year, e.g. from the southern hemisphere) and fitting `c`.
  Planting the seed at a *fixed time* (the season start) as well as a fixed size means a wave's timing
  is then explained by `S0` alone, so "an earlier / steeper season = a more susceptible population"
  holds. Consequence: the *absolute* value of `c` is not interpretable (it is anchored to the
  arbitrary `I0`), whereas *relative* variation in `c` across seasons may still carry meaning
  (care-seeking, testing intensity, severity-driven detection). *Alternative:* fit `I0` per season —
  rejected (confounds the curve's scale and timing with `c` and `S0`).

- **Interpret S0 as RELATIVE across seasons, not absolute.** The absolute level is conditional on the
  fixed R0 and seed (e.g. a lower R0 forces a higher S0); only the ranking / spacing across seasons
  is robust — and that is the quantity of interest.

- **Per season fit S0 and the reporting fraction `c`; per country share the baseline `b` and
  overdispersion `phi`** (plus, for the EKF, the process noise `q_I`). S0 (a rate) and c (a level)
  control different features of the curve, so they separate cleanly; reporting is allowed to vary by
  season (severity / testing / behaviour). *Alternative:* a single shared c — over-constrains size
  vs timing.

- **Off-season activity is an additive observation baseline `b`, not part of the dynamics.** The
  off-season ILI+ floor is non-SIR (sporadic / cross-reacting detections); absorbing it in `b` frees
  the SIR to fit the wave. *Alternative:* restrict the likelihood to the epidemic window — rejected
  as circular (the window depends on the fit we are trying to make).

- **`onset_week` is read off the fitted curve with a simple threshold, and kept as a model feature.**
  For the EKF the curve is the filtered mean, which tracks data and can give early onsets; this is
  accepted as a property of that method rather than something to "fix".

## Method framework

- **A swappable multi-method framework** (`sir_core.R` + `methods/` + `methods_registry.R`): every
  method fits the same per-season panel and returns one common summary table — the input the
  downstream correlation analysis consumes. Adding a method = one new file + one registry line.

- **Primary-analysis method choice is deliberately deferred.** The project is still maturing, so the
  entries below record the *rationale and observed properties* of each method, not a final selection
  of which method(s) lead the primary analysis. All methods are kept and compared; the
  phenomenological vs mechanistic comparison is itself part of the inquiry.

- **Properties of the deterministic SIR (the rigid end of the spectrum).** The deterministic method
  treats the SIR as the truth with residuals as pure noise, which makes it the most rigid: it cannot
  bend to seasons whose shape departs from a single SIR (e.g. the early, sharp FR 2022/23 wave fails
  outright). The **EKF** sits at the other end, admitting genuine process noise (a fitted `q_I` that
  lets the latent state wander off the SIR). Which of these — if either — is carried into the primary
  analysis is left open (see the deferral note); at minimum the deterministic fit is a transparent
  reference / sanity check.

- **EKF process noise is fitted but regularised small, with a tight initial covariance.** This keeps
  S0 identified by the rise rate rather than absorbed by filter freedom. Observed trade-off: the EKF
  draws more between-season S0 contrast (0.70–0.90) than the deterministic (0.73–0.76); the
  descriptive method's rise-rate S0 agrees with the deterministic, indicating the extra EKF spread is
  filter freedom, not raw signal — so read the EKF *ranking*, not the absolute gaps.

- **Take the descriptive (phenomenological) characterisation as a first-class, low-assumption lens —
  and do not over-trust the mechanistic SIR.** A single-population SIR (deterministic or filtered)
  assumes the national ILI+ curve is generated by *one well-mixed epidemic*. National surveillance
  curves are not: they are the **overlay of many local epidemics** that ignite at different times and
  travel across a country, and the aggregate signal is further blurred by reporting delays that
  themselves vary in space and time. Influenza spreads hierarchically, as travelling waves structured
  by human mobility, taking weeks to months to diffuse nationally (Viboud et al. 2006; Gog et al.
  2014). A mechanistic SIR fitted to the aggregate therefore **mis-specifies the data-generating
  process**, and its fitted parameters (e.g. `S0`) can be distorted by aggregation — inviting spurious
  mechanistic interpretation. Following the precautionary, parsimonious principle, we judge it more
  honest to **compare the observable shapes** of the epidemics directly (onset, peak, intensity,
  steepness, burden) and look for predictable patterns than to assume a transmission mechanism we know
  holds only approximately. For an EU agency working closely with the member-state experts who produce
  these data, this is also the more **responsible** stance: phenomenological summaries make minimal,
  transparent assumptions about each country's data instead of imposing a model that may not hold
  uniformly across settings. We therefore treat the phenomenological characterisation as a first-class
  approach alongside the mechanistic fits, and keep all methods; which method(s) lead the primary
  analysis is left open at this stage (see the deferral note), and comparing the phenomenological and
  mechanistic readings is itself part of the inquiry.
  *Method specifics:* smooth each curve and extract AUC, peak height, onset and steepness; steepness
  maps to an implied `S0` via the rise-rate relation, so it sits on the same axis as the SIR methods.
  AUC and peak height are reporting-scale dependent (comparable across seasons *within* a country, not
  across countries); steepness and onset are scale-free.

- **Smoothing = centered moving average, window 4 (not loess).** On regular, single-wave weekly data
  a moving average gives essentially the same features as loess (steepness ~0.31 vs 0.29), so the
  simpler, more transparent option wins. The even "2×4" window (half-weights at the ends → no timing
  shift) preserves the peak better than window 5 while staying smooth. *Alternative:* loess (fine but
  less parsimonious) or a GCV spline (numerically fragile on these short series).

## Data & infrastructure

- **A committed slim panel (`data/slim_flu_iliplus.csv`), loadable in base R.** The susceptibility
  fits run offline from this file (a contiguous weekly grid, seeded from the season start), so they
  need no tidyverse pipeline — fast, dependency-light, reproducible. Countries: DK, FR, IE, HU (each
  with 5 clean seasons; chosen for completeness), set via `params$susc_countries`.

- **Retired the single-season, free-R0 EKF tracking fit.** Superseded by the EKF susceptibility
  method; removing it keeps the repo focused. The shared `.ekf_filter` engine is kept in `sir_core.R`.

- **NA handling.** The moving average is NA-aware and bridges short gaps; a week whose whole window
  is missing stays NA and is skipped by the feature helpers. As more countries bring longer internal
  gaps, the planned step is linear interpolation of internal gaps before smoothing (hook noted in
  `.smooth_curve`); leading / trailing off-season NA can stay (≈ 0). Not built yet — the current
  panel has no problematic gaps.

- **Environment.** CRAN is blocked in the managed environment; dependencies install from Posit
  Package Manager binaries (the setup script sets the repo URL + HTTP user agent) and are restored
  via `renv`.

## References

Methodological and supporting literature for the decisions above. Driver references (dominant
subtype, vaccine coverage, antigenic match / effectiveness, prior immunity) are in `PROJECT_SCOPE.md`
(refs 4–8).

- Viboud C, Bjørnstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT. Synchrony, waves, and spatial
  hierarchies in the spread of influenza. *Science.* 2006;312(5772):447–451.
- Gog JR, Ballesteros S, Viboud C, et al. Spatial transmission of 2009 pandemic influenza in the US.
  *PLoS Computational Biology.* 2014;10(6):e1003635.
- Bjørnstad ON, Finkenstädt BF, Grenfell BT. Dynamics of measles epidemics: estimating scaling of
  transmission rates using a time series SIR model. *Ecological Monographs.* 2002;72(2):169–184.
- Finkenstädt BF, Grenfell BT. Time series modelling of childhood diseases: a dynamical systems
  approach. *Journal of the Royal Statistical Society: Series C (Applied Statistics).*
  2000;49(2):187–205.
- Transmissibility and severity of influenza virus by subtype. *Infection, Genetics and Evolution.*
  2018. https://www.sciencedirect.com/science/article/abs/pii/S1567134818306051

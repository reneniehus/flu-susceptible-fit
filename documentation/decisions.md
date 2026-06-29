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

- **Per season fit S0 and the reporting fraction `c`; share the baseline `b`, overdispersion `phi`
  and (EKF) process noise `q_I` across a country's seasons.** S0 (a rate) and c (a level) control
  different features of the curve, so they separate cleanly, and both are allowed to vary season to
  season (susceptibility; reporting via severity / testing / behaviour). `b`, `phi` and `q_I` are
  instead **shared across the country's seasons**: they are treated as stable country / reporting /
  process properties, and sharing them keeps the parameter count down and **stabilises the fit**.
  *Alternative:* a single shared c (over-constrains size vs timing), or per-season `b` / `phi` / `q_I`
  (more parameters competing with S0 and c).

- **Off-season activity is an additive observation baseline `b`, not part of the dynamics.** The
  off-season ILI+ floor is non-SIR (sporadic / cross-reacting detections); absorbing it in `b` frees
  the SIR to fit the wave. *Alternative:* restrict the likelihood to the epidemic window — rejected
  as circular (the window depends on the fit we are trying to make).

- **`onset_week` is read off the fitted curve with a simple threshold, and kept as a model feature.**
  For the EKF the curve is the filtered mean, which tracks data and can give early onsets; this is
  accepted as a property of that method rather than something to "fix".

- **Neg-binomial-like observation noise** (`Var = mu + mu^2/phi`). Inherited from the Stan generative
  model — but note that Stan fits absolute *counts* whereas we fit a *rate* (ILI+), so a count-style
  variance may not be ideal here. Alternative noise models (e.g. multiplicative / log-scale) should be
  tested — see the planned sensitivity analyses below.

- **Fix the recovery rate `gamma` from a 3-day mean infectious period** (`susc_infectious_period_days`).
  A defensible literature value, fixed for simplicity; estimates span ~2–4 days, and `gamma`
  propagates into `beta = R0*gamma` and the rise-rate scaling, so the value moves the fitted numbers.
  Kept fixed, but flagged for the planned sensitivity analyses below.

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

- **EKF process noise is fitted but regularised small, with a tight initial covariance.** The aim is
  to keep S0 identified by the rise rate rather than absorbed by filter "tracking": a small initial
  covariance `p0` and small process noise `q_I` keep the latent-state covariance small relative to the
  observation noise, so the filter trusts the dynamics over the data. **Caveat (finite sample):** with
  only ~50 weekly observations per season the state covariance still accumulates enough over the
  season that S0 is not fully pinned — seen as the EKF drawing more between-season S0 contrast
  (0.70–0.90) than the deterministic fit (0.73–0.76), while the descriptive rise-rate agrees with the
  deterministic. So read the EKF *ranking*, not the absolute gaps. **This limited-data caveat must be
  stated clearly in the final reporting**, and — if the EKF stays in focus — probed by a sensitivity
  analysis on `p0` / `q_I` (see the planned sensitivity analyses).

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
  *Method specifics:* smooth each curve and extract AUC, peak height, onset and steepness — all
  observed-shape descriptors. None is a clean susceptibility: the observed national rise reflects not
  only transmissibility but also how a country's local epidemics overlay (staggered starts, travel)
  and how reporting is delayed in space and time. AUC, peak height and steepness are therefore
  compared across seasons *within* a country, not across countries (onset timing likewise needs care
  across countries). The method deliberately does **not** map steepness onto an SIR susceptibility —
  that mechanistic-vs-phenomenological mapping is a separate question, out of scope here.

- **Smoothing = centered moving average, window 4 (not loess).** On regular, single-wave weekly data
  a moving average gives essentially the same features as loess (steepness ~0.31 vs 0.29), so the
  simpler, more transparent option wins. The even "2×4" window (half-weights at the ends → no timing
  shift) preserves the peak better than window 5 while staying smooth. *Alternative:* loess (fine but
  less parsimonious) or a GCV spline (numerically fragile on these short series).

## Data & infrastructure

- **Four seasons impacted by the acute COVID-19 pandemic phase are excluded:** 2019/2020, 2020/2021,
  2021/2022 and 2022/2023. They are treated as disrupted by the pandemic and the response to it (NPIs,
  collapsed/atypical influenza circulation, and changes in care-seeking and testing). The analysis
  therefore spans 2014/2015–2018/2019 and 2023/2024–2025/2026 (8 seasons).

- **The two ILI+ sources span different eras and must be combined for a long record.** ERVISS sentinel
  ILI+ (reconstructed as ILI consultation rate × influenza positivity) only reaches back to 2020/21,
  whereas RespiCompass ILI+ runs 2014/15–2023/24. Covering the pre-pandemic seasons (2017/18, 2018/19)
  *and* the most recent ones (2024/25, 2025/26) therefore requires both sources, with **2023/24 as the
  overlap** for cross-checking. Because the two construct ILI+ differently, mixing them along the time
  axis risks a spurious "source/era effect" in the very inter-seasonal comparison of interest; the
  overlap season is used to quantify and adjust for that difference. (See `code/03_report/
  data_availability.R` for the coverage picture.)

- **ERVISS ILI+ is reconstructed RespiCompass-style.** ILI+ = ILI consultation rate × influenza test
  positivity, using **sentinel** positivity (re-derived as detections/tests), **except non-sentinel
  positivity for Malta, Iceland, Croatia, Romania, Latvia and Finland**, which lacked an adequate
  number of weeks of sentinel test positivity — following RespiCompass. This assumes non-sentinel test
  positivity reflects influenza positivity in a similar way as sentinel data would. **Consultation-rate
  units also differ by country:** ERVISS reports per 100 000 population, except Cyprus, Luxembourg and
  Malta (per 100 consultations) and Finland (per 100 000 consultations); to match RespiCompass the
  per-100-consultations countries (CY, LU, MT) are scaled ×1000 onto the per-100 000 basis.

- **The 2023/24 overlap validates the reconstruction, explains the offsets, and anchors the stitch.**
  On 2023/24, **15 of 24 countries match RespiCompass exactly** (cor ≈ 1, ratio ≈ 1). The deviations
  are understood, not method errors:
  - **Consultation-rate units (CY, LU, MT):** off by ×1000 (per-100-consultations vs per-100 000) — a
    deterministic, principled correction; after it LU matches exactly.
  - **Positivity-construction differences (IE ≈ ×0.59, BE ≈ ×1.18, LT):** RespiCompass's influenza
    positivity differs from the re-derived detections/tests by a roughly constant per-country factor;
    once applied the dynamics are identical (IE cor 1.00, BE 0.996; LT 0.91, a few noisy weeks).
  - **LV → ERVISS only (decided).** The RespiCompass LV series is anomalous (flat all winter, a single
    May 2024 spike) while the ERVISS reconstruction is sensible (a January peak) — the source data, not
    the method, is wrong. LV's RespiCompass seasons are therefore dropped and LV is kept on the ERVISS
    reconstruction only (2023/24–2025/26), treated like the single-source countries below.
  - **No clean overlap (NO, ES, SK):** kept on a **single source — whichever has the most non-COVID
    seasons — dismissing the other**: NO → RespiCompass (6 seasons), ES → RespiCompass (5), SK → ERVISS
    (3). No alignment factor is applied to single-source countries (LV likewise: ERVISS only).
  **Stitch:** RespiCompass ≤ 2023/24 + ERVISS reconstruction 2024/25+, with a per-country alignment
  factor from 2023/24 applied to the ERVISS era (= 1 for the 15 exact-match countries; the unit /
  positivity factors above otherwise).

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

## Open questions & planned sensitivity analyses

Choices made for simplicity or by assumption that should be probed downstream, once the analysis
pipeline is in place:

- **Observation-noise model** — test alternatives to the neg-binomial-like variance (the data are
  rates, not counts).
- **Infectious period / `gamma`** — vary around the 3-day default.
- **Fixed `R0` value** — vary around 1.5.
- **Seed `I0` and seed week** — vary the seed size and the season-start week (the latter also tests
  the "early season = more susceptible" reading).
- **EKF `p0` / `q_I` and finite-sample identifiability** — with only ~50 weekly points per season the
  EKF S0 is not fully pinned (see the EKF entry); vary the initial covariance and process-noise
  regularisation to quantify how much of the between-season S0 spread is filter freedom rather than
  signal. To be added at the very end, if the EKF remains in focus.
- **FOLLOW-UP (surveillance/ERVISS experts): why ERVISS and RespiCompass positivity differ for IE,
  BE, LT.** Their ILI+ offset is a *positivity-construction* difference, not a units one — RespiCompass's
  influenza positivity differs from the re-derived ERVISS detections/tests by a roughly constant
  per-country factor (IE ≈ ×0.59, BE ≈ ×1.18; LT matches but is noisy). The exact reason (different
  test denominator, a different/blended positivity source, published vs re-derived) is unknown and
  should be clarified with ERVISS / surveillance colleagues. For now these countries are aligned
  empirically via the 2023/24 overlap factor.

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

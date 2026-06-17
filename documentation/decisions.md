# Design decisions & rationale

The README records *what* the repo does and how to run it; this file records *why* the key
modelling, method and data choices were made, so the reasoning survives beyond commit messages and
chat. Append new decisions as they are made. Keep each entry short: the **decision**, the reason,
and the main alternative considered.

## Inference & modelling

- **Fix R0 from the literature (`susc_R0 = 1.5`); do not fit it.** Jointly fitting R0 and S0 is
  weakly identified (R0 drifts to its prior and S0 follows), so the susceptibility estimate becomes
  meaningless. Fixing R0 lets the wave's *rise rate* `r = gamma*(R0*S0 - 1)` identify S0.
  *Alternative:* fit R0 — rejected (degenerate; R0 drifted to the prior edge, deterministic
  reconstruction poor).

- **Fix the seed `I0` (`susc_seed_i0 = 1e-5`) and plant it at the season start (~Aug).** Encodes a
  constant "flu is always seeded from the southern hemisphere" import. With the seed size and timing
  both fixed, the wave's timing is explained by S0 alone, so "an earlier / steeper season = a more
  susceptible population" holds. *Alternative:* fit I0 per season — rejected (confounds wave timing
  with S0).

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

- **Why we do NOT focus on the deterministic SIR.** The deterministic method (the SIR *is* the truth;
  residuals are pure noise) is rigid: it cannot fit seasons whose shape departs from a single SIR
  (e.g. the early, sharp FR 2022/23 wave fails outright). Real epidemics have process noise, so the
  **EKF** — which admits a fitted `q_I` so the latent state can wander off the SIR — is the preferred
  mechanistic method. The deterministic fit is kept as a transparent *reference / sanity check*, not
  the headline.

- **EKF process noise is fitted but regularised small, with a tight initial covariance.** This keeps
  S0 identified by the rise rate rather than absorbed by filter freedom. Observed trade-off: the EKF
  draws more between-season S0 contrast (0.70–0.90) than the deterministic (0.73–0.76); the
  descriptive method's rise-rate S0 agrees with the deterministic, indicating the extra EKF spread is
  filter freedom, not raw signal — so read the EKF *ranking*, not the absolute gaps.

- **A descriptive (non-mechanistic) method** smooths each curve and extracts AUC, peak height, onset
  and steepness; steepness maps to an implied S0 via the rise-rate relation, so it sits on the same
  axis as the SIR methods and acts as a cross-check. Note: AUC and peak height are reporting-scale
  dependent (comparable across seasons *within* a country, not across countries); steepness and onset
  are scale-free.

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

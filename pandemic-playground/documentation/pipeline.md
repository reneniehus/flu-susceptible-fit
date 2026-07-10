# The simulation pipeline — the maths, stage by stage

Notation: day index `t = 0, 1, …, T−1`; strains `k`; countries `c`; source `X`. `g(·)` is the
discretised generation interval (day-0 mass 0). All stages run under one `set.seed(config$seed)`.

`simulate_pandemic(config)` =
`draw_parameters → simulate_source → simulate_flights → simulate_importations → simulate_local →
observe → assemble`.

## 0. draw_parameters — resolve the config (deterministic)

Turns `epidist` delay objects into day-indexed PMFs by interval censoring,
`w(d) = F(d+½) − F(d−½)` renormalised, with `w(0)=0` for the generation interval (the `cori`
boundary). Collapses an age-structured IFR to an effective scalar `IFR = Σ_a IFR_a weight_a`. Maps
calendar dates (variant introduction, travel ban) to day indices. No random draws here.

## 1. simulate_source — the epidemic at X

A stochastic renewal process with susceptible depletion, per strain `k`:

  Λ_{k,t} = Σ_{s≥1} g(s) I_{k,t−s}                          (total infectiousness of strain k)
  I_{k,t} ~ Poisson( R_{k,t} · (S_{t−1}/N_X) · Λ_{k,t} )  +  seed_{k,t}
  S_t     = S_{t−1} − Σ_k I_{k,t}

`R_{1,t}` is the source's change-point schedule; the variant's `R_{2,t} = R_{1,t}·(1+fitness)`. Seeds:
the day-0 index cases (strain 1) and the variant's introduction on its date (strain 2). Infectious
**prevalence** — the currently-infectious count that drives exports — is incidence convolved with the
survival of the infectiousness profile, `P_t = Σ_s I_{t−s}·(1−CDF_g(s))`; the exportable fraction is
`P_t / N_X`. The variant frequency at source is `I_{2,t}/(I_{1,t}+I_{2,t})`.

## 2. simulate_flights — daily travel volume X → c

  volume_{X→c,t} = baseline_scale · (pop_c/10⁶)^power · connectivity_c · season(t) · ban(t)

`connectivity_c` is a persistent per-country lognormal multiplier; `season(t)=1+amp·cos(2π(t−peak)/365)`;
`ban(t)` multiplies by `travel_ban_factor` after the ban date. Emits the **true** volume matrix and an
**observed** copy blurred by mean-preserving lognormal measurement noise (CV configurable).

## 3. simulate_importations — exports seed the countries

Expected infected arrivals scale with true volume and the source's exportable prevalence:

  E[imports_{c,t}] = volume_{X→c,t} · (P_{X,t}/N_X),     imports_{c,t} ~ Poisson(E[imports_{c,t}])

Each country-day's imports are split into strains by a binomial draw on X's current variant prevalence
fraction, so the variant travels abroad in proportion to connectivity × source dynamics. The **true**
imports (by country, day, strain) are recorded; their *detection* happens later in `observe`.

## 4. simulate_local — the country epidemics

Each country runs the **same** renewal engine as the source, in its own population `N_c`, with its own
change-point `R_{c,t}` (variant `= R_{c,t}·(1+fitness)`), **seeded each day by the imported infections**
from Stage 3 (added as exogenous infections on top of endogenous transmission). Countries differ in
when they intervene, how many imports they receive (connectivity × source prevalence) and their
population — heterogeneous waves from one mechanism.

## 5. observe — the observation model (truth → observed)

Per location, from latent infections `I_t` (total and by strain):

- **Onsets:** `onsets_t = Σ_d I_{t−d}·incubation(d)`.
- **Reported cases (eventual):** `C_t ~ Poisson(onsets_t · ρ_case,t)`, `ρ_case,t` the time-varying
  ascertainment. **Reporting triangle:** each onset day's `C_t` is scattered across report days,
  `triangle[t,·] ~ Multinomial(C_t, onset→report)`; the report day is `t+d`, cells beyond the horizon are
  never seen. `cases_by_report` = anti-diagonal sums (complete to `as_of`); `observed cases_by_onset` =
  within-horizon row sums (**right-truncated**); the eventual `C_t` goes to `truth` (nowcast target).
- **Deaths:** `deaths_t ~ Poisson(IFR · Σ_d I_{t−d}·(incubation⊛onset→death)(d) · ρ_death)`, `ρ_death`
  high. Also `deaths_by_onset` for onset-based CFR.
- **Admissions:** `adm_t ~ Poisson(rate · Σ_d I_{t−d}·(incubation⊛onset→admission)(d))`, well-ascertained.
- **Detected imports:** `Binomial(true imports_{c,t}, surveillance_quality_c)` — the importation signal.
- **Variant typing:** a sequenced subsample `Binomial(C_t, seq_prob)` labelled variant by the onset-day
  strain composition — the frequency series.
- **Flight volumes:** the noisy observed volumes from Stage 2.

## 6. assemble — the two worlds

Returns `list(truth, observed, config, par, latent)`:

- **truth**: `infections`, `onsets`, `cases_by_onset` (eventual), `Rt` (nominal) + `Rt_effective`
  (realized instantaneous, the Cori target), `case_ascertainment`, `ifr`, `imports`, `flight_volumes`,
  `variant_freq`, `capacity`, `prevalence`, `source_size` — all tidy long frames.
- **observed**: `cases_by_onset` (truncated), `cases_by_report`, `reporting_triangle`, `deaths`,
  `admissions`, `detected_imports` (+ surveillance covariate), `variant_cases`, `flight_volumes`, `as_of`.

The realized instantaneous reproduction number stored in truth is
`Rt_effective(t) = [Σ_k R_{k,t}·(S_{t−1}/N)·Λ_{k,t}] / [Σ_k Λ_{k,t}]` — expected new infections over
total infectiousness — which is exactly what a case-based Cori estimate targets.

## References

- Fraser C. Estimating individual and household reproduction numbers in an emerging epidemic. *PLoS ONE.*
  2007;2(8):e758. (renewal / infectiousness profile)
- Cori A, et al. *Am J Epidemiol.* 2013;178(9):1505–1512. (discretised renewal, serial-interval
  discretisation)
- Nishiura H, et al. *PLoS ONE.* 2009;4(8):e6852. (onset→death delay adjustment)

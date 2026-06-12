# flu-susceptible-fit

Fitting **flu ILI+** with a susceptible-reconstruction **SIR model**, in two flavours:

- **Bayesian SIR (Stan)** — `stan/SIR_multiseason_age_vax_2.stan`: an age- and vaccination-structured,
  multi-season SIR fit by HMC (rstan), with scenario projections.
- **Extended Kalman Filter SIR (R)** — `code/01_main_supporting/model_kalman_sir.R`: the same
  generative model (SIR → ILI+, neg-binomial-like noise) fit with an EKF likelihood and
  `optim()` (MAP), entirely in base R — fast and dependency-light.

The repo also carries the full data layer it builds on: ECDC **ERVISS** + **RespiCompass**
surveillance (ILI/ARI, typing/positivity), contact matrices, vaccination and demography, turned
into tidy model-ready tables, with a data-quality / dynamics report.

## Quick start

```r
renv::restore()                       # install pinned dependencies (renv.lock)
source("code/00_main.R")              # build the data and model inputs
```

Then fit the EKF-SIR to one country/season of flu ILI+:

```r
source("code/01_main_supporting/model_kalman_sir.R")
s   <- kalman_sir_series(models_in, "DK", "2023/2024")   # weekly flu ILI+
fit <- fit_kalman_sir(s$value, infectious_period_days = 3)
fit$params                            # R0, S0 (initial susceptible), I0, c, b, phi, qI
kalman_sir_trajectory(fit)            # deterministic curve; n_weeks > length(y) projects forward
```

Demo that fits two countries (Denmark, Greece) and saves a figure:

```sh
Rscript code/04_modelling/fit_kalman_sir_demo.R     # -> output/kalman_sir_fit.png
```

## The model in one paragraph

A latent SIR in population proportions (`dS/dt = -beta S I`, `dI/dt = beta S I - gamma I`) drives
weekly ILI+, which is taken proportional to weekly new infections (`E[y] = c · new infections + b`,
with `b` a small off-season baseline). Observation noise is neg-binomial-like (`Var = mu + mu²/phi`).
The EKF linearises the one-week SIR map to propagate a Gaussian state and compute the likelihood;
weakly-informative priors regularise the classic SIR identifiability degeneracies (R0→∞ with S0→0).
Fitted R0 ≈ 1.1–1.8 and initial-susceptible fractions reproduce observed waves well (correlation
~0.8–0.97 across the test countries); a single-season SIR does not capture secondary strain waves,
by design.

## Tests

```r
Rscript run_tests.R                   # offline, from the committed data/ snapshots
```

Checks the data contracts, the canonical tables, and the model: that the EKF-SIR fit converges
with epidemiologically plausible parameters and reproduces the observed flu ILI+ wave
(`tests/testthat/test-kalman-sir.R`).

## Layout

```
code/00_main.R                 build data + model inputs
code/01_main_supporting/       setup, validate, load_data, gen_model_input, eyeballing,
                               model_kalman_sir (EKF-SIR), model_* scaffolds
code/02_settings/              settings_version0.R (params)
code/03_report/                eyeballing_report.Rmd (data-quality / dynamics report)
code/04_modelling/             fit_kalman_sir_demo.R
stan/                          SIR_multiseason_age_vax_2.stan (Bayesian SIR)
data/                          committed ERVISS / RespiCompass snapshots (offline bootstrap)
output/                        cached data lists + figures (gitignored, regenerated)
tests/testthat/                contract + model tests
documentation/                 data_overview.md, quickstart.md
```

## Reproducibility

`renv.lock` pins all dependencies (R 4.3.3); `renv::restore()` reproduces the environment.
Data is public ECDC ERVISS / RespiCompass surveillance data. See `documentation/quickstart.md`
and `documentation/data_overview.md` for more.

## Note on the Stan model

The Stan model's priors are currently commented out and a few generated-quantities lines need a
fix (see the review notes in commit history); re-enable/repair them before production HMC use. The
EKF-SIR includes the equivalent priors as `optim` penalties and is the ready-to-run path.

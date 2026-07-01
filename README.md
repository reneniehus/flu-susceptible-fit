# flu-susceptible-fit

Reading a per-season sense of population **susceptibility** off **flu ILI+** waves with a
susceptible-reconstruction **SIR model**, through a small **swappable-method framework** plus a
reference Bayesian model:

- **Method framework (R)** — `code/01_main_supporting/sir_core.R` (shared SIR engine) +
  `code/01_main_supporting/methods/` (one file per method) + `methods_registry.R` (registry +
  common per-season summary schema). Every method fixes R0, the infectious period and the seed
  `I0` (from settings) and fits, per season, the **susceptibility** `S0` and a **reporting
  fraction** `c`, with a shared baseline `b` and overdispersion `phi`. Methods so far:
  - `deterministic` — the season's wave *is* a single deterministic SIR; observations are
    overdispersed noise around it. `S0` is identified by the wave's rise rate. Transparent; misfit
    shows up as honest residuals.
  - `ekf` — the same model with an Extended Kalman Filter, admitting **process noise** so the
    epidemic can wander off the deterministic SIR (one shared `q_I` per country, fitted but
    regularised small). The two methods agree on wave shape but the EKF draws more contrast between
    seasons' `S0` — comparing them is exactly what the framework is for.
  - `descriptive` — non-mechanistic: **smooths** each curve (centered moving average) and reads off
    features (AUC, peak height, onset week, steepness) directly. It fits no SIR and reports no `S0` —
    the observed-shape features are compared *within* a country, not mapped onto the SIR
    susceptibility axis (that mechanistic-vs-phenomenological mapping is deliberately out of scope).
- **Bayesian SIR (Stan)** — `stan/SIR_multiseason_age_vax_2.stan`: an age- and vaccination-
  structured, multi-season SIR fit by HMC, with scenario projections (reference model).

The repo also carries the full data layer: ECDC **ERVISS** + **RespiCompass** surveillance
(ILI/ARI, typing/positivity), contact matrices, vaccination and demography, turned into tidy
model-ready tables, with a data-quality / dynamics report.

## Quick start

```r
renv::restore()                       # install pinned dependencies (renv.lock)
source("code/00_main.R")              # build the data and model inputs
```

Fit per-season **susceptibility** across a country's seasons (base R; no pipeline needed — uses the
committed slim panel `data/slim_flu_iliplus.csv`):

```r
source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods_registry.R")

sl  <- load_flu_iliplus_slim("DK")                  # one country's per-season weekly ILI+
fit <- run_method("deterministic", sl, params)      # any registered method, same call
summarise_method_fit(fit)                           # tidy: S0, R_eff, c, peak/onset week, cor, ...
```

`fit$params$S0` is the per-season susceptibility — **interpret the ranking/relative spacing across
seasons**, not the absolute level (which is conditional on the fixed R0/seed).

Run every registered method and save a figure per method:

```sh
Rscript code/04_modelling/fit_methods_demo.R        # -> output/fit_<method>.png
```

## The model in one paragraph

A latent SIR in population proportions (`dS/dt = -beta S I`, `dI/dt = beta S I - gamma I`,
`beta = R0·gamma`) drives weekly ILI+, taken proportional to weekly new infections
(`E[y] = c · new infections + b`). Observation noise is neg-binomial-like (`Var = mu + mu²/phi`).
With R0, the infectious period and the seed `I0` fixed, the wave's **rise rate**
`r = gamma·(R0·S0 − 1)` identifies the per-season susceptibility `S0`; `c` carries the reporting
scale and `b` the off-season baseline. The deterministic method trusts the SIR fully (residuals are
noise); the EKF method adds process noise so the latent state can deviate from it. A single-season,
single-strain SIR does not capture secondary-strain or strongly NPI-shaped seasons, by design.

## Tests

```r
Rscript run_tests.R                   # offline, from the committed data/ snapshots
```

Checks the data contracts and canonical tables, that each method converges with plausible
parameters and reproduces the observed waves (`tests/testthat/test-sir-deterministic.R`), and that
every registered method honours the common summary schema (`tests/testthat/test-methods-registry.R`).

## Layout

```
code/00_main.R                 build data + model inputs, source the methods
code/01_main_supporting/       setup, validate, load_data, gen_model_input, eyeballing,
                               sir_core (shared SIR engine + loaders),
                               methods/ (one file per fitting method),
                               methods_registry (registry + per-season summary schema)
code/02_settings/              settings_version0.R (params, incl. susc_* fixed values)
code/03_report/                eyeballing_report.Rmd (data-quality / dynamics report),
                               data_availability.R (coverage + exclusions heatmap)
code/04_modelling/             build_slim_panel.R (assemble the slim panel), fit_methods_demo.R
                               (run every method, plot + summarise), descriptive_overview.R, ekf_overview.R
code/05_analysis/              the driver analysis: prepare_descriptors.R, analyse_patterns.R,
                               hierarchical_models.R, dominant_subtype.R, bayes_subtype.R,
                               bayes_prior_burden.R, plot_patterns.R, plot_vax_scatter.R
stan/                          SIR_multiseason_age_vax_2.stan (Bayesian SIR)
data/                          committed ERVISS / RespiCompass snapshots + slim_flu_iliplus.csv
output/                        cached data lists + figures (gitignored, regenerated)
tests/testthat/                contract + method tests
documentation/                 quickstart, data_overview, decisions, analysis_strategy,
                               findings_descriptors, documentation.Rmd (see table below)
```

## Documentation

Where each kind of information lives:

| File | Holds |
|---|---|
| `README.md` | what the repo does, quick start, the method framework, layout |
| `documentation/quickstart.md` | how to set up and run |
| `documentation/data_overview.md` | what data is present (`data`, `models_in`, indicators) |
| `documentation/documentation.Rmd` | the model maths / science (SIR, inference, contact matrix) |
| `documentation/decisions.md` | **why** — rationale for key modelling / method / data decisions |
| `documentation/analysis_strategy.md` | the driver analysis — strategy, principles, what we've learned, ranked next steps |
| `documentation/findings_descriptors.md` | results of the descriptor / driver analyses |
| `documentation/external_drivers.md` | externally-sourced drivers (subtype, vaccination, climate, VE) — provenance, values, caveats |
| `documentation/reflections.md` | project-relevance thoughts / interpretations / pondering (also inline `[REFLECTION]` tags) |
| `PROJECT_SCOPE.md` | project scope (in / out of scope), research aim, references |
| inline header comments | why each function / file works the way it does |

New design decisions go in `documentation/decisions.md` (append-only), so the reasoning is kept
with the code rather than only in commit messages.

## Reproducibility

`renv.lock` pins all dependencies (R 4.3.3); `renv::restore()` reproduces the environment.
Data is public ECDC ERVISS / RespiCompass surveillance data. The build chain regenerates every
analysis input from the committed snapshots: `code/00_main.R` writes the model-ready inputs to
`output/models_in.rds` (a gitignored cache), from which `code/04_modelling/build_slim_panel.R`
reproduces the committed panel `data/slim_flu_iliplus.csv` (verified byte-identical) and
`code/05_analysis/prepare_descriptors.R` the descriptor tables the analyses consume. See
`documentation/quickstart.md` and `documentation/data_overview.md` for more.

## Note on the Stan model

The Stan model's priors are currently commented out and a few generated-quantities lines need a
fix (see the review notes in commit history); re-enable/repair them before production HMC use. The
R method framework includes the equivalent priors as `optim` penalties and is the ready-to-run path.

# pandemic-playground

A **pandemic simulation playground for epidemiological testing**: a lean, modular R simulator that
generates **ground-truth** epidemic data for one non-EU/EEA **source location X** and all **EU/EEA
countries** linked by simulated air travel, applies a realistic **observation model**, and returns
**both** the latent truth and the observed surveillance data — so an analytical tool box
(delay-adjusted CFR, importation risk, nowcasting, Rt, forecasting, variant selection, scenarios) can
be **scored against a truth it never saw**.

This is the testbed for a broader project — *modelling tools for a new outbreak and global public
health emergency*. The whole point is the strict separation of two things that real surveillance
conflates:

1. the **latent data-generating process** (infections, transmission, importation, deaths) — what
   actually happened; and
2. the **observation model** (ascertainment, reporting delays, detection, sequencing) — the degraded
   view an analyst is handed.

Everything is **config-driven and reproducible** from one `set.seed`. There are **no external data
files and no epidemiological-package dependencies** — the engine is base R; the tool box adds only the
tidyverse. The tools the brief names ({epiparameter}, EpiEstim, {cfr}, EpiNow2, epinowcast,
deSolve/odin) are **emulated leanly and documented as swap-in extension fronts** (see below).

## Quick start

```r
source("setup.R")                      # engine + analysis tool box (base R + tidyverse)
sim <- simulate_pandemic()             # one synthetic pandemic: truth + observed, from default_config()
sim                                    # a one-screen summary

input <- as_analysis_input(sim)        # the analyst-facing data (tidy schema; real data slots in here)
ga <- growth_analysis(input, "X", window = 20:45)   # e.g. Phase 0: growth rate -> R
score_growth(sim, ga)                  # ... scored against the hidden truth
```

The guided tour runs the whole thing and scores every tool:

```sh
Rscript demo/run_playground.R
```

## The design in one paragraph

`simulate_pandemic(config)` runs a pipeline —
`draw_parameters → simulate_source → simulate_flights → simulate_importations → simulate_local →
observe → assemble` — and returns `list(truth, observed, config, ...)`. A stochastic **renewal**
process (`I_t ~ Poisson(R_t · S_t/N · Σ_s g(s) I_{t-s})`, one engine for the source **and** every
country) generates infections; a fitter **variant** can be introduced at X and travels abroad with the
passengers. **Flights** carry the source's infectious **prevalence** into each country as **imports**,
which seed the local renewal. The **observation model** then degrades the truth exactly as real
surveillance does: infections → onsets (incubation) → a time-varying **ascertained** fraction of
**cases**, scattered across report days into a **reporting triangle** (right-truncated); **deaths** via
the IFR and the onset→death delay; **admissions** for the capacity question; **detected imports**
thinned by each country's surveillance quality; a **sequenced** subsample for the variant frequency;
and noisy **observed flight volumes**. Delay distributions are `epidist` objects — a tiny,
{epiparameter}-shaped class that also accepts a real {epiparameter} object via `as_epidist()`.

## The analysis tool box (scored against truth)

Structured by the phases of a response. Each tool consumes the tidy `observed` schema (or real data in
the same shape) and has a `score_*()` companion that grades it against `truth`.

| Phase | Question | Tool (file) | Lean method → package it emulates |
|---|---|---|---|
| **0** before local intro | Is it spreading, how fast? | `phase0_growth_R` | Poisson-GLM growth → R (Euler–Lotka) → *EpiEstim* |
| **0** | How big is the source really? | `phase0_catchment` | exported-case back-calculation → *arithmetic / Poisson* |
| **0** | Who is under-detecting imports? | `phase0_importation_risk` | imports~flights regression + PI (De Salazar) → *glm* |
| **1** early local growth | How deadly is it? | `phase1_cfr` | delay-adjusted CFR (Nishiura) → *{cfr}* |
| **1** | Growing or shrinking now? | `phase1_rt` | Cori renewal Rt (conjugate Gamma) → *EpiEstim* |
| **1** | How many recent cases, really? | `phase1_nowcast` | reporting-triangle truncation correction → *epinowcast* |
| **2** growth to peak | Will we breach capacity? | `phase2_forecast` | renewal forecast of cases → admissions vs capacity → *EpiNow2* |
| **2** | Are controls working? | `phase2_intervention` | interrupted time series on the growth rate |
| **3** sustained / later waves | Next season, does boosting help? | `phase3_scenarios` | SIRS factorial (RK4), relative + ensembled → *deSolve/odin* |
| **3** | Is a variant taking over? | `phase3_variant_selection` | logistic selection coefficient (binomial `glm`; `nnet::multinom` for >2) |

Run `Rscript demo/run_playground.R` to see all ten, each next to its truth.

## Real data instead of the simulation

Every tool consumes tidy data frames in a fixed schema (`as_analysis_input()`), **never** the sim
object — so real surveillance data in the same shape is a drop-in replacement. Delay distributions
come from the literature as `epidist` (or {epiparameter}) objects. See
`documentation/real_data.md` for the schema, column by column, and worked substitution points.

## Tests

```sh
Rscript run_tests.R                    # 109 self-contained tests; they simulate their own data
```

Cover the delay distributions, config validation, the renewal engine's invariants, the truth/observed
contract (reproducibility, right-truncation, triangle reconciliation), the observation model, and that
every analysis tool recovers the truth it targets. Two known limitations are pinned as tests, not
hidden: importation inflates a country's early Cori Rt, and a fixed-R forecast overshoots at the peak.

## Layout

```
setup.R                     source the engine + tool box in dependency order
R/                          the simulation engine (base R, no dependencies)
  epidist.R                 delay distributions + discretisation ({epiparameter}-shaped)
  config.R                  default_config(), validate_config(), config_subset(), EU/EEA geography
  utils.R                   step schedules, delay convolutions, seasonality
  renewal.R                 the shared stochastic renewal engine (source AND countries)
  draw_parameters.R         config -> discretised PMFs, effective IFR, day indices
  simulate_source.R         Stage 1: the source epidemic at X (one/two strains)
  simulate_flights.R        Stage 2: true + observed daily flight volumes X -> c
  simulate_importations.R   Stage 3: exports ~ Pois(volume x prevalence), by strain
  simulate_local.R          Stage 4: per-country renewal, seeded by imports
  observe.R                 Stage 5: the observation model (truth -> observed)
  assemble.R                simulate_pandemic(): the full pipeline + print method
analysis/                   the tool box (base R + tidyverse)
  analysis_common.R         the real-data seam (schema) + scoring against truth
  phase0_*.R phase1_*.R phase2_*.R phase3_*.R
demo/run_playground.R       the guided, scored tour
tests/testthat/             self-contained tests (the simulator is the fixture)
documentation/              decisions, pipeline maths, tool box, real data, reflections
```

## Documentation

| File | Holds |
|---|---|
| `README.md` | what it is, quick start, the design, the tool box, layout |
| `PROJECT_SCOPE.md` | aim, in / out of scope, references |
| `documentation/pipeline.md` | the simulation maths (renewal, importation, observation) |
| `documentation/analysis_toolbox.md` | the phase-by-phase tool map, what each is scored against, extension fronts |
| `documentation/real_data.md` | the analysis schema and how to substitute real data |
| `documentation/decisions.md` | **why** — the design decisions and their alternatives |
| `documentation/reflections.md` | pondering / open considerations (also inline `[REFLECTION]` tags) |
| inline header comments | why each file / function works the way it does |

## On the emulated packages (and reproducibility)

The brief names a modern epi tool stack ({epiparameter}, EpiEstim, {cfr}, EpiNow2, epinowcast,
deSolve/odin, nnet). Those are the right tools in production; here they are **implemented leanly in
base R** for three reasons: the managed environment blocks CRAN (only the tidyverse, `nnet`, `MASS`,
`mgcv` are present); a self-contained, dependency-light engine is more reproducible and portable; and
implementing each method transparently makes its assumptions — and the biases the playground is built
to expose — visible rather than hidden behind a call. Every tool documents the package it emulates and
the **extension front** where the fuller version drops in (e.g. `as_epidist()` already ingests a real
{epiparameter} object). The engine has **no** dependencies at all; `set.seed(config$seed)` makes every
run byte-identical.

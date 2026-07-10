# A pandemic simulation playground for epidemiological testing

## Background and rationale

When a novel pathogen emerges in one part of the world, the analytical questions a public-health agency
must answer arrive in a predictable sequence — is it spreading, how fast, how deadly, who imports it
next, will hospitals cope, are controls working, is a variant taking over, what might next season hold.
A large tool box exists for each question (delay-adjusted CFR, importation-risk regression, nowcasting,
Rt estimation, renewal forecasting, variant-selection models, scenario projection). But on **real**
emerging-outbreak data the *truth* is never known at the time — so a tool's answer cannot be checked
against what actually happened until much later, if ever.

This playground supplies the missing truth. It simulates a full synthetic pandemic — a source epidemic
that seeds EU/EEA countries through air travel — where **every latent quantity is known exactly**, then
hands the analyst only the **degraded surveillance data** a real response would have. Each analytical
tool can then be run on the observed data and **scored against the ground truth**, so we can see which
methods recover which quantities, under which biases, and how early.

It is deliberately a **testing** instrument, not a forecasting product. It exists so the analytical
tool box can be *built, stress-tested and calibrated* against known answers before it is trusted on the
next real emergency.

## Aim

Provide a lean, modular, reproducible simulator that (1) generates ground-truth epidemic data for a
source location and all EU/EEA countries linked by simulated travel, (2) applies a realistic
observation model, and (3) returns both worlds — so that (4) a phase-structured analytical tool box can
be scored against the truth, and (5) real data can replace the simulated data at any interface.

## Design principles

- **Strict truth/observation separation.** The latent data-generating process and the observation
  model are separate stages; the simulator returns `truth` and `observed` side by side. No tool is ever
  allowed to touch `truth` except through a `score_*()` function.
- **One config, one seed, one pipeline.** Every assumption lives in `default_config()`; a run is
  reproducible from the config plus its seed. Experiments are edits to a copy of the config.
- **Lean but extendable at many fronts.** The core is small and base-R. Extension points are documented
  where they belong: delay distributions accept real {epiparameter} objects; the renewal engine is
  multi-strain and could take age structure; each analysis tool names the fuller package it emulates.
- **Honest about limits.** Where a lean method is biased (ascertainment-inflated growth, importation-
  contaminated Rt, fixed-R forecast overshoot, over-tight CIs), the playground *demonstrates* the bias
  by scoring it, and the limitation is documented and, in two cases, pinned as a test.

## In scope

- A stochastic **renewal** data-generating process for a source location X and the EU/EEA country set,
  linked by a simulated **air-travel importation** model, with an optional fitter **variant**.
- A realistic **observation model**: time-varying case ascertainment, reporting delays and the
  right-truncated reporting triangle, deaths (IFR + onset→death), hospital admissions, surveillance-
  quality-driven import detection, sequenced variant frequencies, noisy flight volumes.
- A **phase-structured analytical tool box** (Phases 0–3) with a **scoring** layer against truth.
- A documented **real-data seam** so real surveillance data can be substituted for the simulation.

## Out of scope (extension fronts, not built)

- **Age-structured transmission.** The renewal engine is single-well-mixed-population per location; an
  age-structured IFR is collapsed to an effective rate. Age structure is a documented extension.
- **Partial cross-immunity between strains.** Strains share one susceptible pool (complete
  cross-immunity). Partial immunity / immune escape is an extension front.
- **Within-country spatial structure.** National curves are one well-mixed epidemic, not a spatial
  overlay of local waves (the very effect the parent flu project cautions about).
- **The production tool stack itself.** {epiparameter}, EpiEstim, {cfr}, EpiNow2, epinowcast,
  deSolve/odin are emulated leanly rather than depended on (CRAN is blocked here; and transparency is a
  goal). Each is a named swap-in point.
- **Behaviourally-driven or mechanistic contact models, economic modules, real geography beyond
  population.** Not modelled.

## Relationship to the parent repository

This sub-project lives on a dedicated branch of the `flu-susceptible-fit` repository and is written to
become its **own repository**: it is self-contained under `pandemic-playground/` with its own README,
docs, tests and dependency story. It shares the parent's engineering culture — tidyverse-flavoured,
readable code, decisions recorded in `documentation/decisions.md`, pondering in `reflections.md`,
references kept with the work — but no code: the parent studies *seasonal influenza susceptibility from
real ILI+ data*, whereas this studies *how well analytical tools recover a known truth on synthetic
emerging-outbreak data*.

## Key references

Methods emulated by the tool box; full citations are also in each tool's file header and in
`documentation/analysis_toolbox.md`.

1. Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and software to estimate time-varying
   reproduction numbers during epidemics. *Am J Epidemiol.* 2013;178(9):1505–1512.
2. Wallinga J, Lipsitch M. How generation intervals shape the relationship between growth rates and
   reproductive numbers. *Proc R Soc B.* 2007;274(1609):599–604.
3. Nishiura H, Klinkenberg D, Roberts M, Heesterbeek JAP. Early epidemiological assessment of the
   virulence of emerging infectious diseases in a population. *PLoS ONE.* 2009;4(8):e6852.
4. De Salazar PM, Niehus R, Taylor A, Buckee CO, Lipsitch M. Identifying locations with possible
   undetected imported cases by using importation predictions. *Emerg Infect Dis.* 2020.
5. Fraser C, et al. Pandemic potential of a strain of influenza A (H1N1): early findings. *Science.*
   2009;324:1557–1561.
6. Höhle M, an der Heiden M. Bayesian nowcasting during the STEC O104:H4 outbreak in Germany.
   *Biometrics.* 2014;70(4):993–1002.
7. Bernal JL, Cummins S, Gasparrini A. Interrupted time series regression for the evaluation of public
   health interventions: a tutorial. *Int J Epidemiol.* 2017;46(1):348–355.
8. Keeling MJ, Rohani P. *Modeling Infectious Diseases in Humans and Animals.* Princeton Univ. Press; 2008.
9. Epiverse-TRACE R packages — {epiparameter}, EpiEstim, {cfr}, EpiNow2, epinowcast — the production
   tools this playground emulates and is designed to swap in.

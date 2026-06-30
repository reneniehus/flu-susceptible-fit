# Analysis strategy, reasoning, and roadmap

*The thinking behind the driver analysis: the strategy, the methodological principles we've arrived at,
what we've learned, and the next steps with their logic. Companion to `PROJECT_SCOPE.md` (the why),
`findings_descriptors.md` (the results), `external_drivers.md` (externally-pulled driver data +
provenance), and `decisions.md` (the modelling/data decisions). Living document — update as it evolves.*

## 1. The strategy in one line

Turn each country-season influenza epidemic into a few interpretable **phenomenological descriptors**
(AUC/burden, peak height, peak week, onset week, steepness), then ask **what explains their variation**
— always **within country**, **borrowing strength across the EU/EEA** by partial pooling, and probing
candidate drivers one at a time.

## 2. Methodological principles (hard-won, now standing rules)

- **Phenomenological, not mechanistic, descriptors.** National ILI+ curves are spatial overlays of many
  local epidemics with space/time-varying reporting (Viboud 2006; Gog 2014). A single-population SIR
  mis-specifies that, so we compare observed *shapes* rather than impute mechanism. (See the
  phenomenological-approach decision in `decisions.md`.)
- **Within country only.** Cross-country comparisons are confounded by latitude, health system, care-
  seeking, age structure and surveillance design. Every association is estimated within country
  (predictors group-mean-centred; country random intercepts).
- **Partial pooling to keep power.** Hierarchical models (country random intercepts, and random slopes
  where estimable) let the **shared EU/EEA slope learn from all countries** while country deviations and
  heterogeneity are estimated. This is how we keep EU/EEA-wide power without cross-country confounding.
- **Scale discipline.** AUC and peak height are reporting-scale dependent (differ ~1000x across
  countries) → only ever compared **within** country (log + demeaning). Steepness, onset week and peak
  week are scale-free → also comparable across countries (which we still avoid for inference).
- **The descriptors are ~three orthogonal axes:** **size** (peak ≈ AUC, within-country r ≈ 0.98),
  **timing** (onset → peak week, r ≈ 0.6, and decoupled from size), and **speed** (steepness, which
  shares almost nothing with size). Onset predicts *when* the peak comes, not *how big*.
- **Confounds we always check:** (a) the **source/era confound** — pre-COVID = RespiCompass, post-COVID
  = ERVISS, so any pre/post difference is inseparable from a measurement change; (b) **reporting
  persistence** — burden magnitude is autocorrelated within country for non-epidemiological reasons, so
  a positive prior→current burden slope can be reporting, a negative one is the cleaner signal.
- **"No detectable effect" ≠ "proven zero."** With ~20–25 countries and few seasons, wide intervals are
  often power limits, stated as such.
- **Adversarial honesty.** Claims are stress-tested (the descriptor analysis went through explicit
  verification that retracted several overclaims); we separate data-supported statements from
  interpretation throughout.

## 3. The mechanistic lens for interpretation

Even though we fit phenomenological descriptors, we read them through epidemiology:
- **Growth rate (steepness)** is the take-off *engine*: r ≈ (R0−1)/generation-interval ≈ transmissibility
  × starting susceptibility × seasonal forcing × **spatial synchrony**. At the national level synchrony
  matters a lot, so steepness is best read as a **synchrony/timing** quantity, not a local-R or burden
  proxy.
- **Burden (AUC/peak)** is a saturating **final-size integral** + duration + reporting scale. It shares
  only susceptibility/R0 with steepness → which is why speed and burden are nearly decoupled.
- **Timing (onset/peak week)** reflects importation timing, susceptibility build-up, school terms and
  climate; the peak is substantially *anchored* (mid-winter), with onset shifting the lead-in.

## 4. Driver scorecard so far

| driver | data | result | why |
|---|---|---|---|
| 65+ vaccination coverage | in-repo (pre-COVID only) | **null** within country | 65+ are transmission *receivers*; the all-age curve is set by higher-contact younger groups; coverage barely varies year-to-year |
| dominant subtype | ERVISS typing (post-COVID only) | **real contrasts but season-confounded** | only 3 post-COVID seasons, each one subtype → subtype ≈ season; can't separate from vaccine match / weather of those years |
| prior-season burden (immunity) | in-repo (lagged) | **weak / null** | subtype rotation breaks cross-protection → total prior burden is a poor immunity proxy; antigenic drift restores susceptibility; reporting persistence masks it |

**The recurring lesson:** the within-country, *season-decoupled* signal for influenza at this national
aggregation is **small**. That is itself a finding, and it sharpens what a "real" driver must look like.

## 5. The core challenge (why drivers keep dissolving)

Candidate drivers fall into three buckets:
- **Continental-per-season** (dominant subtype, EU-wide vaccine effectiveness, broad antigenic match):
  vary by season but are ~shared across countries → **confounded with the season/calendar**. With only a
  handful of usable seasons they collapse into a season effect.
- **Country-level-constant** (demography, latitude, contact structure, school calendar): vary across
  countries but barely within → **between-country only**, which our within-country rule excludes (these
  are *moderators*, not drivers).
- **Genuine country × season** variation, *not* tied to the calendar: **the sweet spot**. Examples:
  local winter climate anomalies; a country's *own* immunity history (same-subtype prior burden);
  subtype *novelty* relative to that country's recent past. These can carry within-country,
  season-decoupled signal.

A "real" driver, by our standard, must be: **within-country, season-decoupled, robust to the source
confound, and replicated across seasons.**

## 6. Next steps (ranked, with logic)

1. **Climate — winter temperature / absolute humidity.** The prime remaining sweet-spot driver: it
   varies year-to-year *within* a country and is mechanistically strong (Tamerius 2013). Hypothesis:
   colder/drier winters → earlier, sharper, larger epidemics. *Needs external reanalysis (ERA5 /
   Copernicus) per country-season; the network proxy may block the ECDC-style domains, so the fetch is
   the risk.* Highest-value single new driver to try.
2. **Accrue more single-source ERVISS seasons (structural fix, low effort, high leverage).** As ERVISS
   adds seasons under one source: (i) the source/era confound shrinks; (ii) subtypes **recur** across
   years, so subtype becomes separable from season; (iii) the proper **same-subtype prior-burden**
   immunity test becomes possible. Much of the current confounding is a *data-quantity* problem that
   time alone fixes — so re-running these analyses each new season is itself a plan.
3. **Subtype novelty / turnover.** Define susceptibility by whether the circulating subtype is *new
   relative to the country's recent history* (drift/shift → larger susceptible pool). Novelty is keyed to
   each country's own past, so it **partly escapes** the season confound that sank the plain-subtype
   analysis. Needs a few more seasons of typing.
4. **Continuous subtype shares** instead of a single dominant label — recovers the country×season
   variation in the typing data that the argmax discards.
5. **Antigenic match / vaccine effectiveness.** Central to the scope (refs 6–7) and more informative than
   coverage. External (ECDC / I-MOVE seasonal VE); partly season-confounded like subtype, but with more
   country/age variation.
6. **Demography / contacts as country-level *moderators*.** We hold the contact matrices and population
   pyramids in `data/`. They are between-country (excluded as within-country drivers) but could explain
   the **heterogeneity** we already see (e.g. why the prior→peak-week effect varies so much across
   countries: random-slope SD 0.79) — i.e. model country random slopes *as a function of* demography.

### Methodological upgrades (orthogonal to the drivers)
- **Full Bayesian fits** (`rstanarm`/`brms`) for proper posteriors, especially tightening small-N
  intervals (subtype, vaccination); our Gibbs sampler already gives this for the simple models and is
  lme4-validated.
- **Joint / multivariate descriptor models** (size, timing, speed together) to estimate a driver's effect
  on the *whole shape* at once and respect the descriptors' covariance.
- **Sub-national / spatial-synchrony decomposition** — the deep fix for the spatial-overlay caveat:
  separate "local transmission" from "national synchrony" so steepness can be interpreted mechanistically.
- **Validate the post-COVID earlier-onset signal** as ERVISS accumulates single-source seasons (it is
  currently entangled with the RespiCompass→ERVISS handover).

## 7. Open questions to revisit

- The **post-COVID earlier-onset** shift (~1.5–2.5 wk) is real-but-confounded with source; needs more
  single-source seasons to confirm.
- **Steepness interpretation** — commit to reading it as spatial synchrony/timing rather than local R?
- The **IE/BE/LT positivity discrepancy** (a data-construction difference vs RespiCompass) is flagged in
  `decisions.md` for follow-up with ERVISS/surveillance experts.
- Avoid the **`rise = peak − onset` definitional coupling** in any composite descriptor.

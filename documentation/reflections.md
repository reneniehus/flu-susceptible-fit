# Reflections, interpretations & pondering

*A dedicated home for project-relevance thoughts — the "why does this matter / how should we think about
this" reasoning that is more speculative than the settled decisions in `decisions.md` or the results in
`findings_descriptors.md`. Entries are dated and may be provisional. Companion to `analysis_strategy.md`.*

**Tagging convention.** Project-relevance thoughts live here as dated entries. For a thought that belongs
next to a specific piece of code or another doc, tag it inline instead: `# [REFLECTION] ...` in R code, or
`[REFLECTION] ...` / a `> [!reflection]` blockquote in Markdown — so `grep -rn "\[REFLECTION\]"` finds them all.

---

## 2026-07 — External driver availability: identifiability, confounding, and mechanism

Prompted by the external-driver availability panel (`code/03_report/driver_availability.R`). Three linked ideas.

### 1. Some drivers vary by country, others by season — this governs what is identifiable

- **Country-varying** drivers: **vaccination coverage**, **climate**. They differ across countries within a
  season.
- **Season-varying** drivers: **dominant subtype**, **vaccine effectiveness (VE)**. They are ~constant
  across countries within a season (continental).

This split is decisive for regression identifiability and confounding:
- Two **season-varying** predictors (subtype, VE) are mutually collinear and collinear with any season
  effect. With only a handful of seasons they **cannot be cleanly separated** — a "subtype effect" and a
  "VE effect" compete to explain the same season-level variation. (Concretely pre-COVID: the two A(H3N2)
  seasons are also the two lowest-VE seasons, so subtype and VE are entangled.)
- **Country-varying** predictors (coverage, climate) vary *within* a season, so their effect is identified
  from **country × season** variation — far more information, and separable from season-level confounders.
- **The sweet spot is a driver that varies by BOTH country and season.** The cleanest example we already
  hold is **effective protection = VE × coverage**: VE is season-level, coverage is country-level, so the
  product varies on both axes and is much better identified than VE alone. Local climate anomalies (once
  pulled per country-season) are the other sweet-spot driver.

### 2. Pre- vs post-COVID is a likely confounder — analyse pre-COVID first

There are plausible **systematic differences between the pre- and post-COVID eras** beyond the virus:
human contact patterns, health-seeking behaviour, and surveillance/reporting standards all changed. Any
pre/post contrast is confounded by these. (It compounds the mechanical RespiCompass→ERVISS source change we
already track.) **Strategy:** fit the analysis on **pre-COVID seasons only** first — one clean era — and
only then, carefully, add the post-COVID seasons with an explicit era control, treating any era term as
soaking up an unknown mixture of behaviour + reporting, not a clean effect.

### 3. A mechanistic lens that stays phenomenological

- **Vaccines act in tandem: effective protection ≈ VE × coverage.** Halving VE and doubling coverage roughly
  cancel, so the meaningful quantity is the **product**, not either alone.
- **EU flu vaccines do not materially change transmission** (coverage is low and concentrated in the
  elderly, who are not the main transmitters). So they do not alter the epidemic's *dynamics* / final size —
  their effect is a **direct, ~multiplicative reduction of the observed burden** in the protected fraction,
  not a change routed through the final-size relationship.
  - *Clarification we should keep straight:* ILI+ is **symptomatic, medically-attended** influenza, and
    I-MOVE VE is precisely VE against medically-attended influenza — so VE × coverage **does** predict a
    (modest) reduction of ILI+ **AUC** in the vaccinated elderly fraction. The "no effect on spread" point
    means there is **no herd/dynamical amplification** of that reduction, so it stays a direct multiplicative
    burden subtraction. (A multiplicative reduction of burden is an **additive shift on log(AUC)** — which is
    why modelling log(AUC) is consistent with this mechanism, not an imposition of one.)
- **Climate acts on transmissibility.** Its effect on burden should therefore resemble the **Rt → final-size**
  relationship (non-linear, saturating), *unlike* the vaccine's direct-subtraction effect.
- **Avoid hard-coding a mechanism by using categorical predictors of AUC.** Subtype is already categorical;
  once climate is pulled, **categorise it** too, and let each category's mean AUC capture whatever
  (non-linear) relationship exists — rather than assuming a functional form. This keeps the analysis
  phenomenological while remaining mechanistically literate.

### Implications for the models (testable predictions)

| Driver | Mechanism | Should predict… | Should NOT predict… |
|---|---|---|---|
| VE × coverage (protection) | direct burden subtraction, no dynamics | **AUC / peak height** (↓) | onset / peak **week** (timing) |
| Dominant subtype | intrinsic transmissibility + population susceptibility | all descriptors (burden & timing) | — |
| Climate (once pulled) | transmissibility → final size | **AUC** (non-linear) | (timing effect uncertain) |

So a clean falsifiable check: **protection should reduce AUC/peak but leave onset/peak-week unmoved.** This
is the logic behind `code/05_analysis/bayes_precovid_ve_subtype.R` (see §"pre-COVID VE + subtype" results in
`findings_descriptors.md`).

**Result (2026-07, prediction NOT confirmed — and that is the point).** In the pre-COVID model the expected
burden reduction was weak/null (protection→AUC -0.08, ns), while the *only* significant protection/VE
association was with **earlier onset** (VE -0.65\*, protection -0.38\*) — a *timing* effect vaccines cannot
cause. So it is **confounding, not mechanism**: with subtype and VE both season-level over 5 seasons, and
65+ coverage nearly constant within country, the VE×coverage "sweet-spot" term was still dominated by its
season-level VE component and inherited the subtype/season timing signal rather than isolating burden. The
lesson stands and sharpens: the vaccine-burden mechanism is **not identifiable at this aggregation/era** —
it needs either (a) far more seasons, or (b) real *within-country* protection variation, or (c) the outcome
that vaccines actually act on (severe outcomes / hospitalisation), not all-age ILI+. The subtype contrasts
(B biggest & latest, H3N2 earliest) are the informative, if season-confounded, signal.

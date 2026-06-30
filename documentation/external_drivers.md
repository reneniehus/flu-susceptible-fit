# External driver data: provenance, values, and caveats

*Externally-sourced candidate drivers pulled to fill the seasons/variables missing from the in-repo data,
for the inter-seasonal-variability analysis. Companion to `analysis_strategy.md` (the why) and
`findings_descriptors.md` (the results). Every value here traces to a named source; uncertainties are kept.*

## How this was gathered (and a hard environment constraint)

This data was assembled via web research. **In this managed environment `WebFetch` and direct data-API
access are blocked (HTTP 403 by the egress policy)** — confirmed against NASA POWER, Open-Meteo, PMC, and
Copernicus-style endpoints. Only `WebSearch` (which returns sourced result summaries) was available. So:

- Literature-based drivers (dominant subtype, vaccine effectiveness, published vaccination-coverage
  figures) **were obtainable** and are recorded below with sources.
- **Reanalysis climate data (per-country-season ERA5 temperature + absolute humidity) is NOT obtainable
  here** — all gridded-data endpoints are blocked. A turn-key pull recipe for an egress-enabled run is given.

Each value below was cross-checked where possible; the four pre-COVID dominant subtypes and two VE seasons
were independently verified against primary sources. Data files live in `data/external/`.

---

## 1. Dominant subtype — pre-COVID seasons (now complete for all 8 panel seasons)

File: `data/external/dominant_subtype_by_season.csv`. "Dominant" = plurality of characterised detections
among A(H1N1)pdm09 / A(H3N2) / B (the WHO/ECDC method, matching `code/05_analysis/dominant_subtype.R`).
Per the user, the continental dominant subtype is assumed ~identical across EU/EEA countries within a season.

| Season | Dominant | B lineage | Basis | Key source |
|---|---|---|---|---|
| 2014/15 | **A(H3N2)** | (Yamagata late) | literature | ECDC season summary; Eurosurveillance ES2015.20.4.21023 |
| 2015/16 | **A(H1N1)pdm09** | B/Victoria | literature | ECDC summary; Eurosurveillance ES.2016.21.13.30184 |
| 2016/17 | **A(H3N2)** (~76% of sentinel) | — | literature | ECDC summary; Eurosurveillance/PMC5883452 |
| 2017/18 | **B** (97% Yamagata) | B/Yamagata | literature | Eurosurveillance/PMC5883452; ECDC AER 2017-18 |
| 2018/19 | **A(H1N1)pdm09** (biphasic) | — | literature | ECDC AER 2018-19; Eurosurveillance ES.2019.24.9.1900125 |
| 2023/24 | A(H1N1)pdm09 | — | ERVISS (in-repo) | dominant_subtype.R; VEBIS 2023/24 |
| 2024/25 | B | B/Victoria | ERVISS (in-repo) | dominant_subtype.R; interim VE 2024/25 |
| 2025/26 | A(H3N2) (subclade K) | — | ERVISS (in-repo) | dominant_subtype.R; Eurosurveillance ES.2026.31.7.2600109 |

**Why this matters:** the original subtype analysis (`bayes_subtype.R`) used only the 3 post-COVID seasons,
each a single subtype — subtype was perfectly confounded with season. Across 8 seasons each subtype now
**recurs in both eras** (A(H1N1): 2015/16, 2018/19, 2023/24 · A(H3N2): 2014/15, 2016/17, 2025/26 · B:
2017/18, 2024/25), so subtype is no longer collinear with the calendar or with the RespiCompass→ERVISS
source change. `code/05_analysis/subtype_8season.R` re-runs the within-country contrast across all 8 seasons.

### 8-season result (within country, partial pooling, net of `era`; SD units, 95% CrI)

| Outcome | H3N2 − H1N1 | B − H1N1 | B − H3N2 |
|---|---|---|---|
| AUC (log) | **+0.18*** | **+0.36*** | **+0.17*** |
| peak height (log) | **+0.21*** | **+0.30*** | +0.09 |
| peak week | **−0.69*** | +0.10 | **+0.79*** |
| onset week | −0.22 | **+0.36*** | **+0.59*** |
| steepness | +0.13 | −0.05 | −0.18 |

(\* = 95% CrI excludes 0; Gibbs matched lme4 to two decimals, R-hat ≤ 1.011.)

**Reading:** **B seasons carry the largest burden** (AUC) and **peak/onset latest**; **A(H3N2) peaks and
onsets earliest**; **A(H1N1) smallest burden**; steepness shows no subtype signal. These are much sharper
and more confident than the 3-season estimates.

**Caveats (important):**
- Subtype is a **season-level** label (continental), so a "subtype effect" still pools only 2–3 seasons per
  subtype — residual season-level confounding remains, just far less than the 1-season-per-subtype version.
- Because the stitch assigns **2023/24 to RespiCompass**, the A(H1N1) seasons are almost all "pre"-sourced
  (subtype×era table: A(H1N1) 63 pre / 2 post), so H1N1-involving contrasts are partly entangled with era;
  **B − H3N2 is the cleanest contrast** (both span both eras). Read H1N1 contrasts with that in mind.

---

## 2. Vaccination coverage, 65+ — post-COVID seasons

File: `data/external/vaccination_coverage_65plus_postcovid.csv`. The in-repo 65+ coverage runs 2012/13–2021/22;
this adds **2023/24 and 2024/25** for the countries that publish it. **2025/26 is essentially unpublished**
as of mid-2026 (only France's all-at-risk figure and an England final report surfaced).

Panel countries with usable figures (10): DK ~76–78, IE ~75, PT ~71, ES 66 (2023/24), FR ~54, IT ~53, NL
(60+, of invited) ~54, BE 51, PL 9.6, SK 11.7. EU aggregate 47.1% (2023). Sources: Eurostat `hlth_ps_immu`,
ECDC season survey, and national agencies (SSI, HSE, SpF, RKI, RIVM, Sciensano, Ministero/Ministerio).

**Comparability caveats (kept in the file's `note` column):**
- **Calendar-year vs season** — Eurostat reports by calendar year (`2023` ≈ the 2023/24 season), national
  agencies by season; close but not identical.
- **Age band differs** — Germany and the Netherlands report **60+**, not 65+ (→ lower); the Netherlands
  figure is **% of invited**, not population coverage (overstates). These are flagged NON-COMPARABLE.
- Several panel countries (AT, CZ, EE, FI, GR, HR, HU, IS, LT, LU, LV, MT, NO, RO, SI) had **no per-country
  65+ figure surface** in search (Eurostat `hlth_ps_immu` holds them, but a direct table pull is blocked).

Given vaccination coverage was a **null** within-country driver in the existing analysis (65+ coverage
varies little year-to-year and 65+ are transmission receivers), this stream is recorded for completeness;
a clean drop-in extension needs the full Eurostat `hlth_ps_immu` table (egress-enabled pull).

---

## 3. Winter climate — constrained to continental in this environment

File: `data/external/climate_winter_continental.csv`. **Per-country-season temperature + absolute humidity is
not obtainable here** (reanalysis APIs blocked). What is recorded is the **continental** European winter
characterisation per season from Copernicus C3S (DJF anomaly vs 1991–2020 where published, plus the regional
pattern): e.g. 2023/24 **+1.44 °C** and 2024/25 **+1.46 °C** (very mild, 2nd/joint-2nd warmest), 2025/26
**+0.09 °C** (near-normal, record-cold Jan 2026); pre-COVID seasons mostly qualitative (2015/16 very mild &
stormy; 2016/17 SE-Europe cold spell; 2017/18 late "Beast from the East"; 2018/19 mild).

**Analytical caveat:** a continental-per-season anomaly is a **season-level** covariate (like subtype) — it
cannot feed the within-country design, which needs **country × season** variation. The signal that makes
climate a "sweet-spot" driver requires the gridded pull below.

### Recommended pull path (egress-enabled environment)
Pull **ERA5 monthly means** — 2-m temperature + 2-m dewpoint + total precipitation + surface pressure —
for all 8 seasons, via the **Copernicus CDS API** (`reanalysis-era5-single-levels-monthly-means`, free CDS
token) or, to avoid the CDS queue, **ARCO-ERA5 / ERA5-on-AWS** Zarr (public buckets, read with
`xarray`+`zarr`). Aggregate to the **Nov–Mar** (or DJF) window, run **zonal statistics over EU/EEA country
polygons** (NUTS/GADM), and **derive absolute/specific humidity from T + dewpoint** (Clausius–Clapeyron).
Cross-check against **E-OBS** (`insitu-gridded-observations-europe`). This yields the cold/dry signal
motivated by Tamerius et al. 2013 and Lowen & Steel — the highest-value remaining within-country driver.

---

## 4. Vaccine effectiveness (VE) — I-MOVE / VEBIS, point estimates + 95% CIs

File: `data/external/vaccine_effectiveness.csv` (41 rows; point estimates and 95% CIs kept verbatim, each
row carrying its source URL). Primary source: the **I-MOVE / I-MOVE+ / VEBIS** multicentre case-control
studies published in **Eurosurveillance**, reporting pooled European VE by (sub)type, age group, and setting
(primary care / hospital). Two seasons independently verified (2014/15 I-MOVE; 2023/24 VEBIS).

Highlights (end-of-season, all ages, primary care, point % [95% CI]):
- **2014/15** — A(H1N1)pdm09 54.2 [31.2, 69.6]; A(H3N2) **14.4 [−6.3, 31.0]** (drifted); B 48.0 [28.9, 61.9]
- **2015/16** — A(H1N1)pdm09 32.9 [15.5, 46.7]
- **2018/19** — A(H3N2) **−1 [−24, 18]** (very low; 46 [8, 68] in 15–64y)
- **2023/24** — A(H1N1)pdm09 52 [44, 59]; A(H3N2) 35 [20, 48]; **B/Victoria 83 [65, 94]**
- **2024/25** — A(H1N1)pdm09 30 [19, 40]; A(H3N2) 38 [26, 49]
- **2025/26** — interim only (mid-2026): influenza A 25–45% (nine studies); France 65+ 28 [17, 37]

Setting and timing matter: hospital ≥65 VE is recorded where available (e.g. 2016/17 A(H3N2) 17 [1, 31]);
cross-study ranges (no single CI) are marked `value_type=study_range`; interim vs end-of-season is flagged.

**Analytical fit:** VE is largely a **season(×subtype)** quantity (some country/age variation), so like
subtype it is mostly a season-level driver — but it is *informative about subtype severity*: e.g. the very
low A(H3N2) VE in 2014/15 and 2018/19 plausibly contributes to those being notable H3N2 seasons, and pairs
naturally with the subtype contrasts in §1. It is more mechanistically meaningful than coverage (a scope ref).

---

## Summary & next steps

| Stream | Coverage pulled | Granularity | Analytical value here |
|---|---|---|---|
| Dominant subtype | **all 8 seasons** (5 pre-COVID added) | continental/season | **High — re-ran the contrast, sharper & de-confounded (§1)** |
| 65+ coverage | 2023/24–24/25, ~10 panel countries | country×season (partial) | Low (null driver; partial; comparability caveats) |
| Climate | 8 seasons | continental only (APIs blocked) | Low here; **high once ERA5 pulled per-country** |
| Vaccine effectiveness | 8 seasons (+CIs) | season(×subtype/age) | Medium — informs subtype severity |

**Highest-value next step:** the **ERA5 per-country-season climate pull** (recipe in §3) in an egress-enabled
environment — the one remaining within-country, season-decoupled "sweet-spot" driver. Second: a full Eurostat
`hlth_ps_immu` table pull to complete 65+ coverage across all 25 panel countries.

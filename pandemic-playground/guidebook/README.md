# The visual guidebook — *Reading an Outbreak*

A staged, data-driven field guide to the analysis tool box: for each phase of an outbreak it shows the
characteristic surveillance picture, the questions that phase poses, the tools that answer them, and
the visual + numerical insight each produces — **scored against the truth the analyst never sees**.

`guidebook.html` is a **self-contained** HTML page (no external requests): the data is embedded, the
charts are hand-built inline SVG, and it renders in the viewer's light or dark theme. Every figure and
number comes from **one seeded default run** of the simulator, so it is fully reproducible.

## Rebuild it

```sh
Rscript guidebook/make_guidebook_data.R   # run one default pandemic; extract each stage's data -> guidebook_data.json
Rscript guidebook/build.R                 # splice the data into the template -> guidebook.html
```

Then open `guidebook/guidebook.html` in a browser (or publish it as an artifact).

## What's here

| File | Role |
|---|---|
| `make_guidebook_data.R` | runs `simulate_pandemic()` + the tool box, captures each stage's surveillance series and tool outputs scored vs truth, writes `guidebook_data.json` |
| `guidebook_template.html` | the page: layout, both-theme design tokens, and a small inline SVG charting engine, with a `__DATA__` placeholder |
| `build.R` | splices the JSON into the template |
| `guidebook_data.json` | the committed data snapshot behind the built page |
| `guidebook.html` | the built, self-contained guidebook |

## The stages it walks

- **Phase 0 — the signal is abroad.** Source growth → R, catchment back-calculation, importation-risk flagging.
- **Phase 1 — it's here, and climbing.** Nowcasting the truncated tail, renewal R<sub>t</sub>, delay-adjusted CFR.
- **Phase 2 — will the hospitals hold?** Depletion-aware capacity forecast, interrupted-time-series on control.
- **Phase 3 — a new variant, and next season.** Logistic selection coefficient, next-season SIRS scenarios.

The data generator depends only on the playground engine + tool box (base R + tidyverse); the build
step is base R. Charting is vanilla JS/SVG in the template — no chart library, no CDN.
